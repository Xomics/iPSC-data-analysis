library(utils)

# Use the current directory as the source
current_directory <- getwd()

# Create "Zipped" folder inside current directory for the archives
destination_directory <- file.path(current_directory, "Zipped")

if (!dir.exists(destination_directory)) {
  dir.create(destination_directory)
}

# Get only top-level folders, excluding "Zipped"
folders <- list.dirs(current_directory, recursive = FALSE, full.names = FALSE)
folders <- folders[folders != "Zipped"]

for (folder_name in folders) {
  folder_path <- file.path(current_directory, folder_name)
  zip_path <- file.path(destination_directory, paste0(folder_name, ".zip"))
  
  # Remove existing zip if it exists
  if (file.exists(zip_path)) {
    file.remove(zip_path)
  }
  
  # Check if folder has contents
  contents <- list.files(folder_path, recursive = TRUE, all.files = TRUE)
  if (length(contents) > 0) {
    zip(zip_path, files = folder_path)
    cat("Zipped:", folder_name, "->", zip_path, "\n")
  } else {
    cat("Skipped empty folder:", folder_name, "\n")
  }
}