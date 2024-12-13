from flask import Flask, request, jsonify, render_template, url_for, send_from_directory, redirect
import requests
from bs4 import BeautifulSoup
import xmltodict
from Bio import SeqIO 
import pytest
import xml.etree.ElementTree as ET
from io import StringIO
from unittest.mock import patch
import subprocess
from flask_sqlalchemy import SQLAlchemy
import os
from datetime import datetime
from collections import Counter
import time
import matplotlib.pyplot as plt
import io
import base64
import re
from xml.etree import ElementTree as ET
from PIL import Image
from pyzbar.pyzbar import decode
from urllib.parse import quote
from urllib.parse import quote_plus
from flask_cors import CORS
import pandas as pd  
from io import BytesIO
from Bio.Restriction import RestrictionBatch, Analysis, AllEnzymes
from Bio.Seq import Seq

app = Flask(__name__)
CORS(app)

BOLD_API_URL = "http://v3.boldsystems.org/index.php/Ids_xml"

blast_endpoint = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
api_key = "4650e1877eb97708a2e62f181bc9749c5508"

NCBI_BLAST_API_URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
IUCN_API_URL = "https://apiv4.iucnredlist.org/api/v4"
IUCN_API_KEY = "UFKmBzriWAiqpaWmt9uo2e8VSovgfTpxKzqx" 

NCBI_TAXONOMY_API_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_SUMMARY_API_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///users.db'
app.config['SECRET_KEY'] = 'supersecretkey'

# Initialize counters
view_count = 0
download_count = 0

UPLOAD_FOLDER = 'uploads/'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


files = []



RESTRICTION_ENZYMES = {
    "EcoRI": {"sequence": "GAATTC", "overhang": "sticky"},
    "BamHI": {"sequence": "GGATCC", "overhang": "sticky"},
    "HindIII": {"sequence": "AAGCTT", "overhang": "sticky"},
    "SmaI": {"sequence": "CCCGGG", "overhang": "blunt"},
    "NotI": {"sequence": "GCGGCCGC", "overhang": "sticky"},
    "PstI": {"sequence": "CTGCAG", "overhang": "sticky"},
    "XhoI": {"sequence": "CTCGAG", "overhang": "sticky"},
    "SpeI": {"sequence": "ACTAGT", "overhang": "sticky"},
}


@app.route("/", methods=["GET"])
def index():
    return render_template("index.html")

template_routes = {
    "login": "login.html",
    'genedna': 'database/genedna.html',
    'proteindata': 'database/protein.html',
    'submission': 'submission.html',
    'dnafulldata1': 'database/dnafulldata1.html',
    'dnafulldata2': 'database/dnafulldata2.html',
    'dnafulldata3': 'database/dnafulldata3.html',
    'dnafulldata4': 'database/dnafulldata4.html',
    'dnafulldata5': 'database/dnafulldata5.html',
    'rnadata': 'database/rna.html',
    'publications': 'publications/publications.html',
    'journal': 'publications/journal.html',
    'books': 'publications/books.html',
    'book1': 'publications/book1.html',
    'articles': 'publications/articles.html',
     "identifyspecies": "analysis/identifyspecies.html",
    "speciessearch": "analysis/blast.html",
    "restriction": "analysis/restriction.html",
    "arlequin":"analysis/arlequin.html",
    "mega":"analysis/mega.html",
    "network1":"analysis/network1.html",
    "dnasp":"analysis/dnasp.html",
    "alignment":"tools/alignment.html",
    "genebarcoding": "tools/barcode.html",
    "graph": "tools/graph.html",
    "visual":"tools/visual.html",
    "collab": "about/collab.html",
    "people": "about/people.html",
    "advisoryboard": "about/advisoryboard.html",
    "history": "about/history.html",
    "about": "about/about.html"

}


@app.route('/process_sequenceblast', methods=['GET', 'POST'])
def process_sequenceblast():
    if request.method == 'GET':
        return render_template('analysis/blast.html')  # Ensure this template exists

    if request.method == 'POST':
        dna_sequence = request.form.get('dna_sequence', '').replace("\n", "").replace("\r", "").strip()
        if not dna_sequence:
            return "Please provide a valid FASTA sequence."


        params = {
            'CMD': 'Put',
            'PROGRAM': 'blastn',
            'DATABASE': 'nt',
            'QUERY': dna_sequence,
            'FORMAT_TYPE': 'XML'
        }

        try:
            response = requests.post(NCBI_BLAST_API_URL, data=params)
            response.raise_for_status()
        except requests.RequestException as e:
            return f"Error communicating with NCBI BLAST API: {e}"

        rid = extract_rid(response.text)
        if not rid:
            return "Failed to retrieve Request ID (RID) from NCBI."

        for _ in range(30):
            try:
                check_params = {
                    'CMD': 'Get',
                    'RID': rid,
                    'FORMAT_TYPE': 'XML'
                }
                status_response = requests.get(NCBI_BLAST_API_URL, params=check_params)
                status_response.raise_for_status()
                if 'Status=WAITING' not in status_response.text:
                    break
            except requests.RequestException as e:
                return f"Error retrieving BLAST results: {e}"
            time.sleep(10)

        if 'Status=WAITING' in status_response.text:
            return "BLAST search is taking too long. Try again later."

      
        results = parse_blast_results(status_response.text)
        if not results:
            return "No results found."

        return render_template('blast_results.html', results=results)

def extract_rid(response_text):
    match = re.search(r'RID\s*=\s*(\w+)', response_text)
    return match.group(1) if match else None

def parse_blast_results(xml_response):
    results = []
    try:
        root = ET.fromstring(xml_response)
        for hit in root.findall(".//Hit"):
            organism_name = hit.findtext("Hit_def")
            accession = hit.findtext("Hit_accession")
            accession_length = hit.findtext("Hit_len")
            scientific_name = extract_scientific_name(organism_name)
            hsp_identity = hit.find(".//Hsp_identity")
            hsp_align_len = hit.find(".//Hsp_align-len")
            hsp_max_score = hit.find(".//Hsp_bit-score")
            hsp_total_score = hit.find(".//Hsp_score")
            hsp_evalue = hit.find(".//Hsp_evalue")
            hsp_query_cover = hit.find(".//Hsp_query-from")

            
            percentage_identity = None
            if hsp_identity is not None and hsp_align_len is not None:
                identity_value = int(hsp_identity.text)
                align_length = int(hsp_align_len.text)
                percentage_identity = (identity_value / align_length) * 100

         
            results.append({
                'scientific_name': scientific_name,
                'organism': organism_name,
                'accession': accession,
                'accession_length': accession_length,
                'identity': f"{percentage_identity:.2f}" if percentage_identity is not None else "N/A",
                'max_score': hsp_max_score.text if hsp_max_score is not None else "N/A",
                'total_score': hsp_total_score.text if hsp_total_score is not None else "N/A",
                'cover_query': hsp_query_cover.text if hsp_query_cover is not None else "N/A",
                'e_value': hsp_evalue.text if hsp_evalue is not None else "N/A"
            })
    except ET.ParseError:
        return []
    return results

def extract_scientific_name(hit_def):
    """
    Extract the scientific name (genus and species) from the `Hit_def` field.
    The scientific name is assumed to be the first two words in the description.
    """

    words = hit_def.split()
    

    if len(words) >= 2:
        return f"{words[0]} {words[1]}"  
    
 
    return hit_def.strip()


@app.route('/query_bold', methods=['POST'])
def query_bold():
    dna_sequence = request.json.get('dna_sequence')
    database = request.json.get('database', 'COX1')

    if not dna_sequence:
        return jsonify({'error': 'No DNA sequence provided!'}), 400

    url = 'http://v3.boldsystems.org/index.php/Ids_xml'
    params = {
        'db': database,
        'sequence': dna_sequence
    }

    try:
        response = requests.get(url, params=params)
        response.raise_for_status()  
    except requests.exceptions.RequestException as e:
        return jsonify({'error': str(e)}), 500

    results = parse_bold_response(response.text)
    return jsonify({'data': results})


def parse_bold_response(xml_data):
    root = ET.fromstring(xml_data)
    results = []
    for match in root.findall('.//match'):
        results.append({
            'Query ID': 'submitted',
            'PID [BIN]': match.find('processid').text + ' [' + (match.find('bin_uri').text or '') + ']',
            'Phylum': match.find('phylum').text,
            'Class': match.find('class').text,
            'Order': match.find('order').text,
            'Family': match.find('family').text,
            'Subfamily': match.find('subfamily').text or '',
            'Genus': match.find('genus').text,
            'Species': match.find('species').text,
            'Indels': match.find('indels').text or '0',
            'ID%': match.find('similarity').text
        })
    return results

def blast_search(query, database, program):
    params = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": query
    }
    response = requests.post(blast_endpoint, params=params)
    return response.text



@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400

    file = request.files['file']

    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400

    if file:
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], file.filename))
        upload_date = datetime.now().strftime('%Y-%m-%d')
        files.append({'filename': file.filename, 'date': upload_date})
        return redirect(url_for('admin_panel'))


@app.route('/delete/<filename>', methods=['POST'])
def delete_file(filename):
    global files
    files = [f for f in files if f['filename'] != filename]

   
    file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    if os.path.exists(file_path):
        os.remove(file_path)

    return redirect(url_for('admin_panel'))


def query_bold_id_engine(sequence):
    api_url = f"http://v3.boldsystems.org/index.php/Ids_xml?db=COX1_SPECIES_PUBLIC&sequence={sequence}"
    response = requests.get(api_url)
    
    if response.status_code == 200:
        return parse_bold_response(response.text)
    else:
        return None


def parse_bold_response(xml_data):
    import xml.etree.ElementTree as ET

    results = []
    root = ET.fromstring(xml_data)
    
    for match in root.findall(".//match"):
        result = {
            'species': match.find('identification').text if match.find('identification') is not None else "Unknown",
            'similarity': match.find('similarity').text if match.find('similarity') is not None else "N/A",
            'specimen': {
                'publicpage': match.find('specimen_url').text if match.find('specimen_url') is not None else "N/A"
            }
        }
        results.append(result)
    
    return {'identifications': {'match': results}}


@app.route('/identify', methods=['POST'])
def identify_sequence():
    dna_sequence = request.form.get('sequence')  
    print("DNA Sequence received:", dna_sequence)  
    
    if not dna_sequence:
        return jsonify({"error": "No DNA sequence provided"}), 400
    
    result = query_bold_id_engine(dna_sequence)
    
    if result:
        print("BOLD API Result:", result)  
        

        matches = result.get('identifications', {}).get('match', [])
        extracted_data = [
            {
                "species": match['species'],
                "similarity": match['similarity'],
                "specimen_page": match['specimen']['publicpage']
            }
            for match in matches
        ]
        
        return jsonify(extracted_data)
    else:
        print("Failed to retrieve data from BOLD API")  # Log failure for debugging
        return jsonify({"error": "Failed to retrieve data from BOLD API"}), 500



@app.route("/<page_name>")
def render_page(page_name):
    if page_name in template_routes:
        return render_template(template_routes[page_name])
    else:
        return "Page not found", 404


@app.route('/view_pdf/<filename>')
def view_pdf(filename):
    global view_count
    view_count += 1  
    return render_template('publications/book1.html', pdf_file=url_for('static', filename=f'pdfs/{filename}'), 
                           view_count=view_count, download_count=download_count)


@app.route('/download_pdf/<filename>')
def download_pdf(filename):
    global download_count
    download_count += 1  
    pdf_dir = os.path.join(app.root_path, 'static/pdfs')  
    return send_from_directory(pdf_dir, filename, as_attachment=True)


def analyze_with_biopython(user_sequence):
    """
    Analyze a DNA sequence with Biopython's restriction enzyme library.
    Returns enzyme cut sites and relevant data.
    """
    dna_sequence = Seq(user_sequence)
    enzymes = RestrictionBatch(AllEnzymes)  # Load all known restriction enzymes
    analysis = Analysis(enzymes, dna_sequence)  # Perform analysis
    cut_sites = []

    for enzyme, sites in analysis.full().items():
        if sites:  # Only include enzymes that actually cut
            ovhg = enzyme.ovhg if enzyme.ovhg is not None else 0  # Handle None values
            if ovhg > 0:
                overhang_type = "sticky"
            elif ovhg == 0:
                overhang_type = "blunt"
            else:
                overhang_type = "sticky (5'/3')"
            
            cut_sites.append({
                "enzyme_name": enzyme.__name__,
                "site_length": len(enzyme.site),
                "sequence": enzyme.site,
                "overhang": overhang_type,
                "cut_positions": sites,
                "frequency": len(sites),
            })
    return cut_sites


def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))

def find_cut_sites(dna_sequence, enzyme_data):
    cut_sites = []
    reverse_dna = reverse_complement(dna_sequence)
    for enzyme, details in enzyme_data.items():
        pattern = details["sequence"]
        positions = [m.start() + 1 for m in re.finditer(pattern, dna_sequence)]
        positions += [m.start() + 1 for m in re.finditer(pattern, reverse_dna)]
        if positions:
            cut_sites.append({
                "enzyme_name": enzyme,
                "site_length": len(pattern),
                "sequence": pattern,
                "overhang": details["overhang"],
                "cut_positions": sorted(positions),
                "frequency": len(positions),
            })
    return cut_sites
  
@app.route('/analyze', methods=['POST'])
def analyze_sequence():
    """
    Analyze the user-provided DNA sequence with restriction enzymes.
    """
    dna_sequence = request.form.get('dna_sequence', '').upper()
    if not re.match("^[ACGT]*$", dna_sequence):
        return jsonify({"error": "Invalid DNA sequence. Only A, C, G, and T are allowed."}), 400
    
    try:
        results = analyze_with_biopython(dna_sequence)
        if not results:
            return jsonify({"message": "No enzymes cut this sequence."})
        return jsonify({"cut_sites": results})
    except Exception as e:
        return jsonify({"error": f"Error analyzing sequence: {str(e)}"}), 500


@app.route('/submit_barcode', methods=['POST'])
def submit_barcode():
    if 'dna_barcode' not in request.files:
        return jsonify({"success": False, "error": "No barcode file uploaded"}), 400

    barcode_file = request.files['dna_barcode']

    try:

        decoded_data = decode_qr_code(barcode_file)

        species, gene_info, similarity_score, accession, database = match_barcode_with_database(decoded_data)

        if species:

            results = {
                "species": species,
                "gene_info": gene_info,
                "similarity_score": similarity_score,
                "accession": accession,
                "database": database
            }
            return jsonify({"success": True, "results": results})
        else:
            return jsonify({"success": False, "error": "No matching result found."})
    except ValueError as e:
        return jsonify({"success": False, "error": str(e)})

@app.route('/download/excel', methods=['GET'])
def download_excel():
 
    df = pd.DataFrame(results)

    output = BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        df.to_excel(writer, index=False, sheet_name='Search Results')
    output.seek(0)

    return send_file(
        output,
        as_attachment=True,
        download_name="Search_Results.xlsx",
        mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
    )

@app.route('/download/csv', methods=['GET'])
def download_csv():

    df = pd.DataFrame(results)
    
   
    output = BytesIO()
    df.to_csv(output, index=False)
    output.seek(0)

    return send_file(
        output,
        as_attachment=True,
        download_name="Search_Results.csv",
        mimetype='text/csv'
    )

if __name__ == "__main__":
    if not os.path.exists(UPLOAD_FOLDER):
        os.makedirs(UPLOAD_FOLDER)
    app.run(debug=True)