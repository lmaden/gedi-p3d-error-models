# save as batch_newlines_to_commas.py
from pathlib import Path

# Set your input directory (use a raw string on Windows to avoid backslash escapes)
INPUT_DIR = Path(r"C:\Users\levim\OneDrive\Desktop\work\gedi\dissertation\chpt1\data\CATIDs_VRIC-141315")

def convert_file(in_path: Path) -> int:
    """Read newline-separated IDs and write a single-line, comma-separated file.
    Returns the number of non-empty lines written."""
    out_path = in_path.with_name(f"{in_path.stem}_comma{in_path.suffix}")

    # Read with 'utf-8-sig' to safely handle a BOM if present
    with in_path.open("r", encoding="utf-8-sig") as fin:
        ids = [line.strip() for line in fin if line.strip()]

    with out_path.open("w", encoding="utf-8") as fout:
        fout.write(",".join(ids))  # no spaces; add ", ".join(...) if you want spaces later

    return len(ids)

def main():
    if not INPUT_DIR.exists():
        raise SystemExit(f"Input directory not found: {INPUT_DIR}")

    txt_files = sorted(p for p in INPUT_DIR.glob("*.txt") if not p.stem.endswith("_comma"))
    if not txt_files:
        print("No .txt files to process (or only *_comma.txt files present).")
        return

    total_files = 0
    for p in txt_files:
        try:
            n = convert_file(p)
            print(f"Converted {p.name} -> {p.stem}_comma.txt  ({n} IDs)")
            total_files += 1
        except Exception as e:
            print(f"!! Skipped {p.name} due to error: {e}")

    print(f"Done. Processed {total_files} file(s).")

if __name__ == "__main__":
    main()
