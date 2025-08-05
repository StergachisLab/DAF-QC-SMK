import panel as pn
import base64
from pathlib import Path

# This script was created with the assistence of claude (claude.ai)

pn.extension(template='material')

sample_name = snakemake.params.sample_name
file_paths = snakemake.input.pdfs
regions = snakemake.params.regions
output_file = snakemake.output.dashboard

#sample_name = "Sample Name"
#file_paths= ["/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test_samp/qc/reads/plots/htt_test_samp.chr3_179228176_179236561.bias.reads.pdf", "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test_samp/qc/reads/plots/htt_test_samp.chr4_3073138_3075853.duplication_groups.pdf"]
#regions= ["chr3:179228176-179236561", "chr4:3073138-3075853"]
#output_file= "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test_samp/qc/testdash2.html"


def parse_filenames(filenames, regions):
    pdfs = {}

    plot_type = {
        "duplication_groups": "Duplicate Distribution",
        "strandtype": "Strand Type",
        "bias": "DddA Sequence Bias",
        "deam_rate": "Overall Deamination Rates",
        "mut_rate": "Mutation Rates",
        "targeting_plot": "Targeting Metrics",
        "other": "Other"
        }
    
    
    for pdf_path in filenames:
        pdf_path = Path(pdf_path)
        if not pdf_path.exists() or pdf_path.suffix.lower() != ".pdf":
            continue

        with open(pdf_path, 'rb') as f:
            pdf_data = f.read()
            
        pdf_base64 = base64.b64encode(pdf_data).decode()

        name = pdf_path.stem

        region = next((r for r in regions if r in name), "All")

        category = next((cat for cat in plot_type.keys() if cat in name), "other")

        readtype = next((rt for rt in ["reads", "consensus"] if rt in name), "")


        plot_name = f"{readtype.capitalize()}: {plot_type[category]}" if readtype else plot_type[category]

        # Load associated text file, if it exists

        text_content = ""
        text_path = pdf_path.with_suffix('.txt')
        if text_path.exists():
            try:
                with open(text_path, 'r', encoding='utf-8') as f:
                    text_content = f.read()
            except Exception as e:
                text_content = f"Error reading text file: {e}"


        pdfs[pdf_path.stem] = {
            'base64': pdf_base64,
            'region': region,
            'display_name': plot_name,
            'path': str(pdf_path),
            'text_content': text_content
        }

    return pdfs

def find_text_file(pdf_path):
    """Find text file that matches the PDF file"""
    pdf_stem = pdf_path.stem
    
    for text_file in text_files:
        text_path = Path(text_file)
        if text_path.stem == pdf_stem:
            return text_path
    return None

def create_standalone_html(pdfs, sample_name, output_file):
    """Create a standalone HTML file with embedded PDFs"""
    
    # Group PDFs by region
    regions = {}
    for key, pdf in pdfs.items():
        reg = pdf['region']
        if reg not in regions:
            regions[reg] = []
        regions[reg].append((key, pdf))
    
    # Get first PDF for initial display
    first_pdf_key = list(pdfs.keys())[0] if pdfs else None

    # Create the complete HTML
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>PDF Dashboard - {sample_name}</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        body {{
            margin: 0;
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: #f5f5f5;
        }}
        
        .container {{
            display: flex;
            height: 100vh;
        }}
        
        .sidebar {{
            width: 240px;
            background: white;
            border-right: 1px solid #ddd;
            padding: 20px;
            overflow-y: auto;
            box-shadow: 2px 0 5px rgba(0,0,0,0.1);
        }}
        
        .sidebar h1 {{
            color: #2c3e50;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        
        .sidebar h3 {{
            color: #34495e;
            margin: 20px 0 10px 0;
            font-size: 18px;
        }}
        
        .sidebar h4 {{
            color: #7f8c8d;
            margin: 15px 0 8px 0;
            font-size: 16px;
            border-bottom: 1px solid #ecf0f1;
            padding-bottom: 5px;
        }}
        
        .sidebar ul {{
            list-style: none;
            padding: 0;
            margin: 0 0 15px 0;
        }}
        
        .sidebar li {{
            margin: 5px 0;
        }}
        
        .sidebar a {{
            display: block;
            padding: 8px 12px;
            text-decoration: none;
            color: #2c3e50;
            border-radius: 5px;
            transition: background 0.3s;
            cursor: pointer;
        }}
        
        .sidebar a:hover {{
            background: #ecf0f1;
        }}
        
        .sidebar a.active {{
            background: #3498db;
            color: white;
        }}
        
        .main-content {{
            flex: 1;
            padding: 20px;
            overflow: hidden;
            display: flex;
            gap: 20px;
        }}
        
        .pdf-section {{
            flex: 2;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            display: flex;
            flex-direction: column;
        }}
        
        .text-section {{
            flex: 1;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            display: flex;
            flex-direction: column;
            max-width: 240px;
        }}
        
        .pdf-container {{
            height: calc(100vh - 40px);
            display: flex;
            flex-direction: column;
        }}
        
        .pdf-header {{
            padding: 15px 20px;
            border-bottom: 1px solid #ecf0f1;
            background: #f8f9fa;
            border-radius: 8px 8px 0 0;
        }}
        
        .pdf-header h2 {{
            margin: 0;
            color: #2c3e50;
            font-size: 20px;
        }}
        
        .text-header {{
            padding: 15px 20px;
            border-bottom: 1px solid #ecf0f1;
            background: #f8f9fa;
            border-radius: 8px 8px 0 0;
        }}
        
        .text-header h2 {{
            margin: 0;
            color: #2c3e50;
            font-size: 18px;
        }}
        
        .pdf-viewer {{
            flex: 1;
            padding: 0;
        }}
        
        .pdf-embed {{
            width: 100%;
            height: 100%;
            border: none;
        }}
        
        .text-content {{
            flex: 1;
            padding: 20px;
            overflow-y: auto;
            font-family: 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
            font-size: 12px;
            line-height: 1.4;
            background: #fafafa;
            white-space: pre-wrap;
            word-wrap: break-word;
        }}
        
        .pdf-section-inner {{
            display: none;
            height: 100%;
        }}
        
        .pdf-section-inner.active {{
            display: flex;
            flex-direction: column;
        }}
        
        .text-container {{
            display: none;
            height: 100%;
        }}
        
        .text-container.active {{
            display: flex;
            flex-direction: column;
        }}
        
        .no-selection {{
            display: flex;
            align-items: center;
            justify-content: center;
            height: 100%;
            color: #7f8c8d;
            font-size: 18px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="sidebar">
            <h1>{sample_name}</h1>
            <h3>Regions</h3>"""
    
    # Add sidebar navigation
    for region, pdfs_in_region in sorted(regions.items()):
        html_content += f"<h4>{region}</h4><ul>"
        for key, pdf in pdfs_in_region:
            active_class = "active" if key == first_pdf_key else ""
            html_content += f"""
                <li><a class="pdf-link {active_class}" data-pdf="{key}">{pdf['display_name']}</a></li>
            """
        html_content += "</ul>"
    
    html_content += """
        </div>
        
        <div class="main-content">
            <div class="pdf-section">
                <div class="pdf-container">
    """
    
    # Add PDF sections
    for i, (key, pdf) in enumerate(pdfs.items()):
        active_class = "active" if key == first_pdf_key else ""
        html_content += f"""
                    <div id="{key}" class="pdf-section-inner {active_class}">
                        <div class="pdf-header">
                            <h2>{pdf['region'] + ' ' + pdf['display_name']}</h2>
                        </div>
                        <div class="pdf-viewer">
                            <embed src="data:application/pdf;base64,{pdf['base64']}" 
                                   class="pdf-embed" type="application/pdf">
                        </div>
                    </div>
            """
    
    # Add "no selection" message if no PDFs
    if not pdfs:
        html_content += """
                    <div class="no-selection">
                        <p>No PDFs available</p>
                    </div>
        """
    
    html_content += """
                </div>
            </div>
            
            <div class="text-section">
    """
    
    # Add text containers
    for i, (key, pdf) in enumerate(pdfs.items()):
        active_class = "active" if key == first_pdf_key else ""
        text_content = pdf['text_content'].replace('<', '&lt;').replace('>', '&gt;')
        html_content += f"""
                <div id="text-{key}" class="text-container {active_class}">
                    <div class="text-header">
                        <h2>Details</h2>
                    </div>
                    <div class="text-content">{text_content}</div>
                </div>
        """
    
    # Add default text if no PDFs
    if not pdfs:
        html_content += """
                <div class="no-selection">
                    <p>No text content available</p>
                </div>
        """
    
    html_content += """
            </div>
        </div>
    </div>
    
    <script>
        // Handle PDF navigation
        document.addEventListener('DOMContentLoaded', function() {
            const pdfLinks = document.querySelectorAll('.pdf-link');
            const pdfSections = document.querySelectorAll('.pdf-section-inner');
            const textContainers = document.querySelectorAll('.text-container');
            
            pdfLinks.forEach(link => {
                link.addEventListener('click', function(e) {
                    e.preventDefault();
                    
                    const targetPdf = this.getAttribute('data-pdf');
                    
                    // Remove active class from all links and sections
                    pdfLinks.forEach(l => l.classList.remove('active'));
                    pdfSections.forEach(s => s.classList.remove('active'));
                    textContainers.forEach(c => c.classList.remove('active'));
                    
                    // Add active class to clicked link and corresponding section
                    this.classList.add('active');
                    const targetSection = document.getElementById(targetPdf);
                    const targetTextContainer = document.getElementById('text-' + targetPdf);
                    
                    if (targetSection) {
                        targetSection.classList.add('active');
                    }
                    if (targetTextContainer) {
                        targetTextContainer.classList.add('active');
                    }
                });
            });
        });
    </script>
</body>
</html>"""
    
    # Write HTML to file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    return output_file

# Text files variable
text_files = [
    # Add paths to your text files here
    # Example: "/path/to/text/file1.txt"
]

# Process the data
regions = [reg.replace(":", "_").replace("-", "_") for reg in regions]
pdfs = parse_filenames(file_paths, regions)

if pdfs:
    # Create the standalone HTML file
    output_path = create_standalone_html(pdfs, sample_name, output_file)