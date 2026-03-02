"""
PDF to Images Converter.
Converts specified pages of the Sugar PDF to PNG images for Gemini processing.
"""

import argparse
from pathlib import Path
from pdf2image import convert_from_path

from config_Sugar1974 import PDF_PATH, PAGE_IMAGES_DIR, PDF_DPI


def convert_pdf_to_images(
    pdf_path: Path = PDF_PATH,
    output_dir: Path = PAGE_IMAGES_DIR,
    first_page: int | None = None,
    last_page: int | None = None,
    dpi: int = PDF_DPI,
) -> list[Path]:
    """
    Convert PDF pages to PNG images.

    Args:
        pdf_path: Path to the PDF file.
        output_dir: Directory to save images.
        first_page: First page to convert (1-indexed, inclusive). None = first page.
        last_page: Last page to convert (1-indexed, inclusive). None = last page.
        dpi: Image resolution.

    Returns:
        List of paths to generated image files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"📄 Converting PDF pages to images at {dpi} DPI...")
    print(f"   PDF: {pdf_path.name}")

    kwargs = {"dpi": dpi, "fmt": "png", "thread_count": 4}
    if first_page is not None:
        kwargs["first_page"] = first_page
    if last_page is not None:
        kwargs["last_page"] = last_page

    images = convert_from_path(str(pdf_path), **kwargs)

    start_num = first_page if first_page else 1
    image_paths = []

    for i, image in enumerate(images):
        page_num = start_num + i
        image_path = output_dir / f"page_{page_num:03d}.png"
        image.save(str(image_path), "PNG")
        image_paths.append(image_path)
        print(f"   ✅ Saved page {page_num} → {image_path.name}")

    print(f"   📁 {len(image_paths)} page image(s) saved to {output_dir}/")
    return image_paths


def get_existing_images(
    output_dir: Path = PAGE_IMAGES_DIR,
    first_page: int | None = None,
    last_page: int | None = None,
) -> list[Path]:
    """
    Get already-converted page images from disk.

    Returns:
        Sorted list of image paths, filtered by page range if specified.
    """
    all_images = sorted(output_dir.glob("page_*.png"))

    if first_page is None and last_page is None:
        return all_images

    filtered = []
    for img_path in all_images:
        # Extract page number from filename: page_005.png → 5
        try:
            page_num = int(img_path.stem.split("_")[1])
        except (IndexError, ValueError):
            continue

        if first_page is not None and page_num < first_page:
            continue
        if last_page is not None and page_num > last_page:
            continue
        filtered.append(img_path)

    return filtered


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Sugar PDF pages to images")
    parser.add_argument("--pdf", type=str, default=str(PDF_PATH), help="Path to PDF")
    parser.add_argument("--pages", type=int, nargs=2, metavar=("START", "END"),
                        help="Page range (1-indexed, inclusive)")
    parser.add_argument("--dpi", type=int, default=PDF_DPI, help=f"Image DPI (default: {PDF_DPI})")
    args = parser.parse_args()

    first_page = args.pages[0] if args.pages else None
    last_page = args.pages[1] if args.pages else None

    convert_pdf_to_images(
        pdf_path=Path(args.pdf),
        first_page=first_page,
        last_page=last_page,
        dpi=args.dpi,
    )
