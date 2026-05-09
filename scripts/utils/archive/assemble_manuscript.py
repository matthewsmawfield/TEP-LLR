import os
import re
from pathlib import Path

def strip_html_tags(text):
    """Remove HTML tags and clean up whitespace."""
    # Remove script and style elements
    text = re.sub(r'<(script|style)[^>]*>[^<]*</(script|style)>', '', text, flags=re.IGNORECASE)
    # Remove HTML comments
    text = re.sub(r'<!--.*?-->', '', text, flags=re.DOTALL)
    # Replace common tags with markdown equivalents
    text = re.sub(r'<h2[^>]*>\s*', '\n## ', text, flags=re.IGNORECASE)
    text = re.sub(r'<h3[^>]*>\s*', '\n### ', text, flags=re.IGNORECASE)
    text = re.sub(r'<h4[^>]*>\s*', '\n#### ', text, flags=re.IGNORECASE)
    text = re.sub(r'<h5[^>]*>\s*', '\n##### ', text, flags=re.IGNORECASE)
    text = re.sub(r'</h[2-5][^>]*>\s*', '\n\n', text, flags=re.IGNORECASE)
    text = re.sub(r'<p[^>]*>\s*', '\n', text, flags=re.IGNORECASE)
    text = re.sub(r'</p>\s*', '\n\n', text, flags=re.IGNORECASE)
    # Lists
    text = re.sub(r'<ul[^>]*>\s*', '\n', text, flags=re.IGNORECASE)
    text = re.sub(r'</ul>\s*', '\n', text, flags=re.IGNORECASE)
    text = re.sub(r'<ol[^>]*>\s*', '\n', text, flags=re.IGNORECASE)
    text = re.sub(r'</ol>\s*', '\n', text, flags=re.IGNORECASE)
    text = re.sub(r'<li[^>]*>\s*', '- ', text, flags=re.IGNORECASE)
    text = re.sub(r'</li>\s*', '\n', text, flags=re.IGNORECASE)
    # Other inline elements
    text = re.sub(r'<strong[^>]*>\s*', '**', text, flags=re.IGNORECASE)
    text = re.sub(r'</strong>\s*', '**', text, flags=re.IGNORECASE)
    text = re.sub(r'<em[^>]*>\s*', '*', text, flags=re.IGNORECASE)
    text = re.sub(r'</em>\s*', '*', text, flags=re.IGNORECASE)
    text = re.sub(r'<sup[^>]*>\s*', '^', text, flags=re.IGNORECASE)
    text = re.sub(r'</sup>\s*', '^', text, flags=re.IGNORECASE)
    text = re.sub(r'<sub[^>]*>\s*', '~', text, flags=re.IGNORECASE)
    text = re.sub(r'</sub>\s*', '~', text, flags=re.IGNORECASE)
    # Tables
    text = re.sub(r'<table[^>]*>\s*', '\n', text, flags=re.IGNORECASE)
    text = re.sub(r'</table>\s*', '\n', text, flags=re.IGNORECASE)
    text = re.sub(r'<tr[^>]*>\s*', '', text, flags=re.IGNORECASE)
    text = re.sub(r'</tr>\s*', '\n', text, flags=re.IGNORECASE)
    text = re.sub(r'<th[^>]*>\s*', '| **', text, flags=re.IGNORECASE)
    text = re.sub(r'</th>\s*', '** |', text, flags=re.IGNORECASE)
    text = re.sub(r'<td[^>]*>\s*', '| ', text, flags=re.IGNORECASE)
    text = re.sub(r'</td>\s*', ' |', text, flags=re.IGNORECASE)
    # Remove all remaining tags
    text = re.sub(r'<[^>]+>', '', text)
    # Clean up whitespace: remove leading indentation from each line
    lines = text.split('\n')
    cleaned_lines = []
    for line in lines:
        # Strip leading whitespace but preserve empty lines
        cleaned_line = line.lstrip()
        cleaned_lines.append(cleaned_line)
    text = '\n'.join(cleaned_lines)
    # Normalize multiple newlines
    text = re.sub(r'\n{3,}', '\n\n', text)
    text = re.sub(r'^\s+', '', text)
    return text.strip()

def assemble_manuscript():
    project_root = Path(__file__).resolve().parent.parent.parent
    components_dir = project_root / 'site' / 'components'
    output_file = project_root / '17-TEP-LLR-v0.1-Lucknow.md'
    
    # Read version info
    version_file = project_root / 'VERSION.json'
    version_info = "v0.1 (Lucknow)"
    if version_file.exists():
        import json
        try:
            with open(version_file, 'r') as f:
                vdata = json.load(f)
                version_info = vdata.get('version', version_info)
        except (ValueError, KeyError, IOError) as e:
            pass
    
    # Define the order of components
    component_files = [
        "0_abstract.html",
        "1_introduction.html",
        "2_theory.html",
        "3_methodology.html",
        "4_results.html",
        "5_discussion.html",
        "6_conclusion.html",
        "7_references.html",
        "8_reproducibility.html"
    ]
    
    # Build content with proper header like Jakarta paper
    content_parts = []
    content_parts.append("# Temporal Equivalence Principle: Lunar Laser Ranging and the Nordtvedt Effect")
    content_parts.append("")
    content_parts.append("**Matthew Lukin Smawfield**")
    content_parts.append(f"Version: {version_info}")
    content_parts.append("DOI: [10.5281/zenodo.19446029](https://doi.org/10.5281/zenodo.19446029)")
    content_parts.append("")
    content_parts.append("---")
    content_parts.append("")
    
    for filename in component_files:
        filepath = components_dir / filename
        if filepath.exists():
            print(f"Adding {filename}...")
            with open(filepath, 'r') as f:
                html_content = f.read()
                markdown_content = strip_html_tags(html_content)
                content_parts.append(markdown_content)
                content_parts.append("")
        else:
            print(f"Warning: {filename} not found.")
    
    # Join with proper spacing
    content = '\n'.join(content_parts)
    # Final cleanup
    content = re.sub(r'\n{3,}', '\n\n', content)
    
    with open(output_file, 'w') as f:
        f.write(content)
    
    print(f"Successfully assembled manuscript to {output_file}")

if __name__ == "__main__":
    assemble_manuscript()
