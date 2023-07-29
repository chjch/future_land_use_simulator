from pathlib import Path
import zipfile


def create_zip_with_files(
    target_dir, sources, root_dir="future_land_use_simulator"
):
    # Convert the target_dir to a Path object
    target_path = Path(target_dir)
    # Check if the destination folder exists, if not create it
    target_path.mkdir(parents=True, exist_ok=True)

    zip_file_name = target_path / f"{root_dir}.zip"
    # remove the zip file previously built if already exists
    zip_file_name.unlink(missing_ok=True)

    # Create the zip file and add the files from all source folders to it
    with zipfile.ZipFile(zip_file_name, "w") as zip_file:
        for source in sources:
            source_path = Path(__file__).resolve().parent / source
            if source_path.is_file():
                relative_path = (
                    Path(root_dir) /
                    source_path.relative_to(source_path.parent)
                )
                zip_file.write(source_path, str(relative_path))
            elif source_path.is_dir():
                for file_to_zip in source_path.glob("**/*"):
                    if file_to_zip.is_file():
                        relative_path = (
                            Path(root_dir) /
                            file_to_zip.relative_to(source_path.parent)
                        )
                        zip_file.write(file_to_zip, str(relative_path))

    print(f"All files from {', '.join(sources)} added to {zip_file_name}")


if __name__ == "__main__":
    source_dirs = ["gui", "lib", "data"]
    source_files = [
        "__init__.py",
        "metadata.txt",
        "plugin.py",
        "resources.py",
        "resources.qrc",
        "icon.png"
    ]
    create_zip_with_files(
        Path(__file__).resolve().parent / "build", source_dirs + source_files
    )
