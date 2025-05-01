create_template = function(name, type = c("R", "python"), overwrite = FALSE, out_dir = ".") {
  type = match.arg(type)
  src = system.file("templates", paste0(name, ".", type), package = "sequoiar")
  if (src == "") stop("Template not found.")
  dest = file.path(out_dir, paste0(name, ".", type))
  if (file.exists(dest) && !overwrite) stop("File exists. Use overwrite = TRUE.")
  file.copy(src, dest, overwrite = overwrite)
  message("Template created at ", dest)
}