# プロジェクトのルートディレクトリを自動で見つけ、プロジェクト内のパス管理関数を読み込み、ルートを登録する関数
source_edm_paths <- function(start = getwd()) {
  find_root <- function(path) {
    current <- normalizePath(path, winslash = "/", mustWork = TRUE)

    repeat {
      marker <- file.path(current, ".edm_code_root")
      if (file.exists(marker)) {
        return(current)
      }

      parent <- dirname(current)
      if (identical(parent, current)) {
        stop("Could not find edm_code root.", call. = FALSE)
      }
      current <- parent
    }
  }

  root <- find_root(start)
  source(file.path(root, "R", "project_paths.R"))
  set_project_root(root)
  invisible(root)
}
