import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# 1) GeoJSON を読み込む（ファイル名を適宜変更）
gdf = gpd.read_file("japan_prefectures.geojson")  # または .shp / URL

# 2) 都道府県名カラム名を確認（英語/日本語で differ するので要確認）
# print(gdf.columns)
# print(gdf.head())

# 例：都道府県名が 'prefecture' や 'NAME_1' や 'nam_ja' のどれかに入っていることが多いです。
# 以下は 'prefecture' または 'NAME_1' を試す。適宜修正してください。
if 'prefecture' in gdf.columns:
    name_col = 'prefecture'
elif 'NAME_1' in gdf.columns:
    name_col = 'NAME_1'
elif 'name' in gdf.columns:
    name_col = 'name'
else:
    # とりあえず index を使う
    name_col = gdf.columns[0]

# 3) 投影（描画のために Web メルカトルを利用すると見た目が良い）
gdf = gdf.to_crs(epsg=3857)

# 4) 沖縄県だけ取り出す（GeoJSON に日本語で「沖縄県」なら '沖縄県' を使う）
okinawa = gdf[gdf[name_col].astype(str).str.contains("沖縄")].copy()
if okinawa.empty:
    # 別表記（例 'Okinawa'）を試す
    okinawa = gdf[gdf[name_col].astype(str).str.contains("Okinawa", case=False)]
if okinawa.empty:
    raise SystemExit("沖縄県の行が見つかりません。name_col を確認してください。")

# 5) 全図を作成
fig, ax = plt.subplots(figsize=(8, 11))
# 都道府県を薄い色で塗る（カテゴリ色をつけたい場合は column に名前指定）
gdf.plot(ax=ax, linewidth=0.4, edgecolor='k', column=name_col, cmap='Pastel1', legend=False)
# 沖縄を強調
okinawa.plot(ax=ax, color='orange', edgecolor='k')

ax.set_title("Japan (prefectures) — Okinawa highlighted", fontsize=14)
ax.axis('off')

# 6) 沖縄の拡大 inset（全体図の左下などに挿入）
axins = inset_axes(ax, width="40%", height="30%", loc='lower left',
                   bbox_to_anchor=(0.05, 0.05, 0.4, 0.3), bbox_transform=ax.transAxes)

# inset の表示範囲は沖縄の bounds（少し余白を足す）
minx, miny, maxx, maxy = okinawa.total_bounds
pad = 20000  # 表示余白（メートル単位、投影がメートル系のため）
axins.set_xlim(minx - pad, maxx + pad)
axins.set_ylim(miny - pad, maxy + pad)

# 背景に日本全体を薄く描く（座標が合うように同じ crs）
gdf.plot(ax=axins, color='lightgray', edgecolor='none')
okinawa.plot(ax=axins, color='orange', edgecolor='k')

# 7) 離島（多部位ポリゴン）を分解して centroid に注記を付ける
# explode() で multi-part を分離（geopandas >= 0.8）
okinawa_parts = okinawa.explode(index_parts=True, ignore_index=True)
# centroid を計算（注意: centroid は投影済みであることが望ましい）
okinawa_parts['centroid'] = okinawa_parts.geometry.centroid

# ラベル（番号）を配置してから、どの番号がどの島か確認して名前を追加する方法が現実的
for i, row in okinawa_parts.iterrows():
    c = row['centroid']
    axins.annotate(f"{i+1}", xy=(c.x, c.y), xytext=(3,3), textcoords='offset points', fontsize=9, fontweight='bold')

# 右下に小説明（例：番号と島名の対応は手動で追記）
ax.text(0.02, 0.02, "Inset: Okinawa prefecture (numbers → island parts)", transform=ax.transAxes, fontsize=9)

# 8) 保存
plt.savefig("japan_okinawa_inset.png", dpi=300, bbox_inches='tight')
plt.show()

