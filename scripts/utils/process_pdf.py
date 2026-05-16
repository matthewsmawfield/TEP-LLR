#!/usr/bin/env python3
"""Unified PDF Processing Script
Compresses PDF and embeds comprehensive metadata in one operation.

This script processes the TEP-LLR manuscript PDF ("The Temporal Equivalence 
Principle: Lunar Laser Ranging and the Nordtvedt Effect") by compressing it 
for web distribution and embedding complete academic metadata for proper indexing and citation.

Usage:
    python process_pdf.py <input_pdf> [--quality ebook|printer|prepress|default]
    
Example:
    python process_pdf.py site/public/docs/17-TEP-LLR-v0.1-Lucknow.pdf --quality ebook
"""

import subprocess
import sys
import os
from pathlib import Path
import argparse
import tempfile


def compress_pdf(input_path, output_path, quality='ebook'):
    """Compress PDF using Ghostscript."""
    quality_settings = {
        'screen': '/screen',      # 72 dpi
        'ebook': '/ebook',        # 150 dpi
        'printer': '/printer',    # 300 dpi
        'prepress': '/prepress',  # 300 dpi, color preserving
        'default': '/default'
    }

    if quality not in quality_settings:
        raise ValueError(f"Quality must be one of: {', '.join(quality_settings.keys())}")

    gs_quality = quality_settings[quality]

    # Get original size
    original_size = os.path.getsize(input_path)

    # Compress using Ghostscript
    cmd = [
        'gs',
        '-sDEVICE=pdfwrite',
        '-dCompatibilityLevel=1.4',
        f'-dPDFSETTINGS={gs_quality}',
        '-dNOPAUSE',
        '-dQUIET',
        '-dBATCH',
        f'-sOutputFile={output_path}',
        input_path
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True)
        compressed_size = os.path.getsize(output_path)
        reduction = ((original_size - compressed_size) / original_size) * 100

        return {
            'original_mb': original_size / (1024 * 1024),
            'compressed_mb': compressed_size / (1024 * 1024),
            'reduction_pct': reduction
        }
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Ghostscript compression failed: {e.stderr.decode()}")


def embed_metadata(pdf_path, metadata):
    """Embed metadata into PDF using exiftool."""
    cmd = ['exiftool']

    # Add all metadata fields
    for key, value in metadata.items():
        cmd.extend([f'-{key}={value}'])

    # Overwrite original
    cmd.extend(['-overwrite_original', pdf_path])

    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Exiftool metadata embedding failed: {e.stderr.decode()}")


def verify_metadata(pdf_path, expected_fields):
    """Verify metadata was embedded correctly."""
    cmd = ['exiftool'] + [f'-{field}' for field in expected_fields] + [pdf_path]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError:
        return None


def main():
    parser = argparse.ArgumentParser(
        description='Compress PDF and embed metadata in one operation'
    )
    parser.add_argument('input_pdf', help='Path to input PDF file')
    parser.add_argument(
        '--quality',
        choices=['screen', 'ebook', 'printer', 'prepress', 'default'],
        default='ebook',
        help='Compression quality (default: ebook)'
    )
    parser.add_argument(
        '--doi',
        default='10.5281/zenodo.19446029',
        help='DOI to embed in metadata'
    )

    args = parser.parse_args()

    input_path = Path(args.input_pdf).resolve()

    if not input_path.exists():
        print(f"Error: File not found: {input_path}")
        sys.exit(1)

    print(f"Processing: {input_path}")
    print(f"Quality: {args.quality}")
    print()

    # Step 1: Compress PDF
    print("Step 1: Compressing PDF...")
    with tempfile.NamedTemporaryFile(suffix='.pdf', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        stats = compress_pdf(str(input_path), tmp_path, args.quality)

        # Replace original with compressed version
        os.replace(tmp_path, str(input_path))

        print(f"  Original:    {stats['original_mb']:.2f} MB")
        print(f"  Compressed:  {stats['compressed_mb']:.2f} MB")
        print(f"  Reduction:   {stats['reduction_pct']:.1f}%")
        print()

    except Exception as e:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        print(f"Error during compression: {e}")
        sys.exit(1)

    # Step 2: Embed metadata
    print("Step 2: Embedding metadata...")

    # Paper metadata - must match manuscript, CITATION.cff, and zenodo.txt
    metadata = {
        # Core identification
        'Title': 'Temporal Equivalence Principle: Lunar Laser Ranging and the Nordtvedt Effect',
        'Author': 'Matthew Lukin Smawfield',
        'Creator': 'Matthew Lukin Smawfield',

        # Scientific abstract with key results
        'Subject': (
            'The Temporal Equivalence Principle (TEP) is a scalar-tensor theory in which proper time is a '
            'dynamical field φ that couples universally to all matter via a conformal metric. The coupling '
            'strength is density-dependent through a TEP suppression mechanism. Suppression operates via the '
            'continuous spatial profile of the time field (Temporal Topology), in which high ambient density '
            'in deep potential wells suppresses the local field gradient (Temporal Shear). The degree of '
            'gradient suppression scales with gravitational compactness (Φ/c²). TEP preserves the Weak Equivalence '
            'Principle through universal conformal coupling, but predicts violation of the Strong Equivalence '
            'Principle (SEP) via compactness-dependent suppression. This work tests for SEP violation using '
            'Lunar Laser Ranging data, analysing 25,445 cleaned O-C residuals from five international stations spanning '
            '35 years (1984-2019). The headline estimand is precision-weighted full-systematic regression on the full '
            'cleaned sample: η = -3.91 × 10⁻⁴ ± 5.63 × 10⁻⁵ at 6.94σ (6.78σ cluster-robust). Unweighted full-systematic '
            'OLS yields η = -4.06 × 10⁻⁴ ± 6.58 × 10⁻⁵ at 6.17σ as a sensitivity upper bound. Cook\'s-Distance excision on '
            'the unweighted full-systematic row (N = 23,837) returns η = -3.87 × 10⁻⁴ ± 4.95 × 10⁻⁵ at 7.82σ as a '
            'secondary leverage diagnostic. AR(1) GLS on the full-systematic design remains sign-consistent at 4.51σ. '
            'DE430 residuals provide sign-consistent '
            'cross-ephemeris validation at η = -7.03 × 10⁻⁴ ± 2.12 × 10⁻³. '
            'The observed negative sign establishes that gravitational potential suppression dominates in the '
            'Earth-Moon system, providing the first physical explanation for the unexplained ~1 cm synodic '
            'residual documented by Müller & Nordtvedt (1998).'
        ),

        # Keywords for indexing
        'Keywords': (
            'Temporal Equivalence Principle; TEP; Lunar Laser Ranging; LLR; '
            'Equivalence Principle; Nordtvedt Effect; Strong Equivalence Principle; '
            'Post-Newtonian; Scalar-Tensor Gravity; Geometric Suppression; '
            'Gravitational Physics; Earth-Moon System; INPOP19a; DE430'
        ),

        # Production metadata
        'Producer': 'TEP-LLR Research Project - Version 0.1 (Lucknow)',

        # Rights and identifiers
        'Copyright': 'Creative Commons Attribution 4.0 International License (CC BY 4.0)',

        # Dates
        'CreationDate': '2026:05:14 00:00:00',
        'ModifyDate': '2026:05:14 00:00:00',

        # XMP Dublin Core metadata (exiftool uses these prefixes)
        'XMP-dc:Creator': 'Matthew Lukin Smawfield',
        'XMP-dc:Title': 'Temporal Equivalence Principle: Lunar Laser Ranging and the Nordtvedt Effect',
        'XMP-dc:Description': 'TEP test via Lunar Laser Ranging Nordtvedt effect measurement',
        'XMP-dc:Rights': 'CC BY 4.0',
        'XMP-dc:Identifier': f'doi:{args.doi}',
        'XMP-dc:Source': 'https://github.com/matthewsmawfield/TEP-LLR',
        'XMP-dc:Publisher': 'Zenodo',
        'XMP-dc:Date': '2026-05-14',
        'XMP-dc:Type': 'Preprint',
        'XMP-dc:Format': 'application/pdf',
        'XMP-dc:Language': 'en',

        # PRISM (Publishing Requirements for Industry Standard Metadata)
        'XMP-prism:DOI': args.doi,
        'XMP-prism:URL': 'https://github.com/matthewsmawfield/TEP-LLR',
        'XMP-prism:VersionIdentifier': '0.1',
        'XMP-prism:PublicationName': 'TEP Research Series',

        # PDF/A metadata
        'XMP-pdfaid:Part': '1',
        'XMP-pdfaid:Conformance': 'B'
    }

    try:
        embed_metadata(str(input_path), metadata)
        print("  Metadata embedded successfully")
        print()

    except Exception as e:
        print(f"Error during metadata embedding: {e}")
        sys.exit(1)

    # Step 3: Verify
    print("Step 3: Verifying metadata...")
    verification = verify_metadata(
        str(input_path),
        ['Title', 'Author', 'Subject', 'Keywords', 'Creator', 'Copyright']
    )

    if verification:
        print("  ✓ Metadata verified")
        print()
        print("Verification output:")
        print(verification)
    else:
        print("  ⚠ Could not verify metadata")

    print()
    print(f"✓ Processing complete: {input_path}")
    print(f"  Final size: {os.path.getsize(input_path) / (1024 * 1024):.2f} MB")


if __name__ == '__main__':
    main()
