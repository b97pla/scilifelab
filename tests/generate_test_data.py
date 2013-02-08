
import tempfile
import random
import datetime
import string
import os

def generate_fc_barcode():
    """Generate a flowcell barcode on the format ABC123CXX
    """
    return "{}CXX".format("".join([random.choice("{}{}".format(string.ascii_uppercase,string.digits)) for i in xrange(6)]))

def generate_run_id(fc_barcode=generate_fc_barcode()):
    """Generate a run identifier
    """
    return "{}_{}_0{}_{}{}".format(datetime.date.today().strftime("%y%m%d"),
                                           generate_instrument(),
                                           random.randint(101,999),
                                           random.choice("AB"),
                                           fc_barcode)

def generate_project():
    """Generate a project name on the format J.Doe_90_11
    """
    return "{}.{}_{}_{}".format(
        random.choice(string.ascii_uppercase),
        "".join([random.choice(string.ascii_lowercase) for i in xrange(6)]).capitalize(),
        random.randint(90,99), 
        random.randint(11,99))

def generate_sample(project_id=random.randint(500,999)):
    """Generate a sample name on the format P123_123B
    """
    return "P{}_{}{}_index{}".format(
                    project_id,
                    random.randint(101,199),
                    ["","A","B","C","F"][random.randint(0,4)],
                    random.randint(1,24))

def generate_nucleotide_sequence(**kwargs):
    """Generate a random string of specified length and nucleotides sampled from the specified alphabet
    """
    return ''.join([random.choice(kwargs.get('alphabet','ACGT')) for i in xrange(kwargs.get('sequence_length',101))])

def generate_barcode(len=6):
    """Generate a nucleotide barcode
    """
    return generate_nucleotide_sequence(sequence_length=len)

def generate_sample_file(sample_name=generate_sample(), barcode=None, lane=random.randint(1,8), readno=1):
    """Generate a Casava 1.8+-style sample file name
    """
    if barcode is None: 
        barcode = generate_barcode()
    return "{}_{}_L00{}_R{}_001.fastq.gz".format(sample_name,
                                                barcode,
                                                lane,
                                                readno)

def generate_instrument(prefix="SN"):
    """Generate an instrument identifier
    """
    return "{}{}".format(prefix,str(random.randint(101,9999)))

def generate_fastq_header(**kwargs):
    """Generate a fastq header on the CASAVA 1.8+ format:
    @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x- pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>
    """
    
    return "@{} {}".format(':'.join([kwargs.get('instrument',generate_instrument()),
                                      str(kwargs.get('run_number',random.randint(101,9999))),
                                      kwargs.get('fcid',generate_fc_barcode()),
                                      str(kwargs.get('lane',random.randint(1,8))),
                                      str(kwargs.get('tile',random.randint(1101,2316))),
                                      str(kwargs.get('xpos',random.randint(1001,9999))),
                                      str(kwargs.get('ypos',random.randint(1001,9999)))]),
                            ':'.join([str(kwargs.get('read',random.randint(1,2))),
                                      kwargs.get('is_filtered','N'),
                                      str(kwargs.get('control',0)),
                                      kwargs.get('index',generate_barcode())]))

def generate_quality_sequence(**kwargs):
    """Generate a random quality string with the specified length and in the specified span with the specified ascii offset
    """
    return ''.join([chr(random.randint(kwargs.get('qlow',0),kwargs.get('qhigh',40))+kwargs.get('ascii_offset',33)) for i in xrange(kwargs.get('sequence_length',101))])

def generate_samplesheet_data(barcode=generate_fc_barcode(), no_lanes=8, no_projects=2, no_samples=2, genome_build='hg19'):
        
    csv_data = []
    
    # Generate data for each lane
    for lane in range(1,(no_lanes+1)):
        # Create some projects
        for i in range(no_projects):
            project_name = generate_project()
            # Create some samples
            for j in range(no_samples):
                sample_name = generate_sample()
                # Append the generated data to the csv_data
                csv_data.append([barcode,
                                 lane,
                                 sample_name,
                                 genome_build,
                                 generate_barcode(),
                                 project_name,
                                 'C',
                                 'R',
                                 'O',
                                 project_name])
    return csv_data

def generate_run_samplesheet(barcode=generate_fc_barcode(), dst_file=None):
    return _write_samplesheet(generate_samplesheet_data(barcode),dst_file)

def _parse_samplesheet(samplesheet):
    rows = []
    with open(samplesheet) as inh:
        for row in inh:
            rows.append(row.strip().split(","))
    return rows

def _write_samplesheet(csv_data, dst_file=None):
    
    header = ["FCID",
              "Lane",
              "SampleID",
              "SampleRef",
              "Index",
              "Description",
              "Control",
              "Recipe",
              "Operator",
              "SampleProject"]
    
    if dst_file is None:
        fh, dst_file = tempfile.mkstemp(suffix=".csv", prefix="SampleSheet")
        os.close(fh)
    with open(dst_file,"w") as out_handle:
        for row in ([header] + sorted(csv_data, key=lambda x: x[1])):
            out_handle.write(",".join([str(item) for item in row]))
            out_handle.write("\n")
    return dst_file
   
def generate_fastq_record(**kwargs):
    """Generate a fastq record. If the optional parameter 'pair' is True, generate a pair sequence as well
    """
    record = [generate_fastq_header(**kwargs),
              generate_nucleotide_sequence(**kwargs),
              '+',
              generate_quality_sequence(**kwargs)]
    if kwargs.get('pair',False):
        record[0] = record[0].replace(' 2:',' 1:')
        record += [record[0].replace(' 1:',' 2:'),
                   generate_nucleotide_sequence(**kwargs),
                   '+',
                   generate_quality_sequence(**kwargs)]
    return record