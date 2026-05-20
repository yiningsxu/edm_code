find_edm_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    marker <- file.path(current, ".edm_code_root")
    if (file.exists(marker)) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find edm_code root. Start R inside edm_code or a project subdirectory.", call. = FALSE)
    }
    current <- parent
  }
}

edm_root <- function() {
  find_edm_root()
}

workspace_root <- function() {
  override <- Sys.getenv("MICROBIOME_DYNAMICS_ROOT", unset = NA_character_)
  if (!is.na(override) && nzchar(override)) {
    return(normalizePath(override, winslash = "/", mustWork = TRUE))
  }
  normalizePath(file.path(edm_root(), ".."), winslash = "/", mustWork = TRUE)
}

data_path <- function(...) {
  file.path(workspace_root(), "data", ...)
}

result_path <- function(...) {
  file.path(workspace_root(), "result", ...)
}

project_file <- function(...) {
  file.path(edm_root(), ...)
}

set_project_root <- function(root = workspace_root()) {
  setwd(root)
  invisible(root)
}
