import fcntl
import hashlib
import sys
import os

try:
    import pymongo
    from pymongo import MongoClient
except:
    pass
    
TEMPLATE_SEQUENCE_DARPIN = {
    'BCL2L2' : 'DLGKKLLEAARAGQDDEVRILMANGADVNATDASGLTPLHLAATYGHLEIVEVLLKHGADVNAIDIMGSTPLHLAALIGHLEIVEVLLKHGADVNAVDTWGDTPLHLAAIMGHLEIVEVLLKHGADVNAQDKFGKTAFDISIDNGNEDLAEILQK'
}

DEFAULT_SEQUENCE = TEMPLATE_SEQUENCE_DARPIN['BCL2L2']

ALPHABETS = {
    'default': {
        '+': [ 'K', 'R', 'H' ],
        '-': [ 'D', 'E' ],
        'p': [ 'S', 'T', 'N', 'Q' ],
        'u': [ 'A', 'V', 'I', 'L', 'M' ],
        'a': [ 'F', 'Y', 'W' ],
    }
}

REVERSED_ALPHABETS = {}
for ab in ALPHABETS:
    REVERSED_ALPHABETS[ab] = {}
    for k,v in ALPHABETS[ab].items():
        for r in v:
            REVERSED_ALPHABETS[ab][r] = k
        
N_SHORTUID_CHARS = 12
# MAX_MOLECULES = 31600 # make sure these values are in sync with the C-programs
MAX_MOLECULES = 3000 # make sure these values are in sync with the C-programs
MAX_PREDICTIONS = 1000000000

TARGET_PDBS = {
    'BCL2L2' : [ '4k5aA', '4k5bC' ]
}

CONTROL_PDBS = {
    'BCL2L2' : ['4k5aB', '4k5bB' ]
}

CONTROLS = {
    'BCL2L2' : [
        '300025b43772', # 4k5a (10.3 nM)
        '4b1975acce85'  # 4k5b (643 pM)
    ]
}

HITSET_CHAINS = {
}

homedir = os.environ['HOME']
darpinsdir = os.path.join(homedir,'Projects','Darpins')
modelsdir = os.path.join(darpinsdir,'models')
targetsdir = os.path.join(darpinsdir,'targets')
controlsdir = os.path.join(darpinsdir,'controls')
dockingdir = os.path.join(darpinsdir,'docking')
dockinginputdir = os.path.join(dockingdir,'input')
dockingoutputdir = os.path.join(dockingdir,'output')
datadir = os.path.join(darpinsdir,'data')
tmpdir = os.path.join(darpinsdir,'tmp')
workdir = os.path.join(darpinsdir,'work','%d' % os.getpid())


def get_filter_data(args):
    data = {}
    if args.type:
        data.update({'type': args.type })
    if args.parent:
        data.update({'parent': args.parent })
    if data and args.include_parent:
        data = { "$or" : [ data, { 'uid': args.parent } ] }
    return data
        
def get_parent_uid(designs,parent):
    if parent == '': return parent
    col = 'shortuid' if len(parent)==N_SHORTUID_CHARS else 'uid'
    d = designs.find_one({ col: parent })
    if not d: return None
    return d['uid']
    
def is_prepared(mol2file):
    prmtopfile = mol2file.replace('mol2','prmtop')
    prepared = os.path.isfile(mol2file) and os.stat(mol2file).st_size > 0 and \
               os.path.isfile(prmtopfile) and os.stat(prmtopfile).st_size > 0
    if not prepared:
        sys.stderr.write("WARNING: '%s' is not prepared\n" % mol2file)
    return prepared

def hash(s):
    return hashlib.sha224(bytes(s,'utf-8')).hexdigest()
    
def hash_seq(seq):
    return hash(''.join(seq))

def get_mongo_collection(coll,test=False):
    mongodb = 'mongodb://localhost:27017/darpins'
    client = MongoClient(mongodb)
    darpins = client.darpins
    return getattr(darpins,'%s%s' % (coll,('_test' if test else '')))

def get_mongo_predictions(test=False):
    return get_mongo_collection('predictions',test=test)

def get_mongo_designs(test=False):
    return get_mongo_collection('designs',test=test)

def get_designs_from_mongodb(test=False):
    designs = get_mongo_designs(test=test)
    return designs.find().sort([("uid", pymongo.ASCENDING)])

def get_designs_from_dbfile(file):
    designs = []
    lines = open(file).readlines()
    for line in lines:
        tmp = line.rstrip().split('\t')
        design = { 'shortuid': tmp[0], 'seq': tmp[1] }
        designs.append(design)
    return designs

def get_designs(dbfile=None,test=False):
    return get_designs_from_dbfile(dbfile) if dbfile else get_designs_from_mongodb(test=test)

def build_folder(folder):
    if not os.path.isdir(folder):
        os.makedirs(folder)
        
def build_folders():
    build_folder(modelsdir)
    build_folder(tmpdir)
    build_folder(dockingdir)
    build_folder(dockinginputdir)
    build_folder(dockingoutputdir)
    build_folder(workdir)
    build_folder(datadir)

def append_lines_to_file(file,lines):
    with open(file, "a") as out:
        fcntl.flock(out, fcntl.LOCK_EX)
        for line in lines:
            out.write(line)
        out.flush()
        fcntl.flock(out, fcntl.LOCK_UN)

def read_sites_map_from_file(file):
    sites = {}
    if not os.path.isfile(file): return sites
    lines = []
    with open(file, "r") as in_:
        fcntl.flock(in_, fcntl.LOCK_EX)
        lines = in_.readlines()
        fcntl.flock(in_, fcntl.LOCK_UN)
    for line in lines:
        try:
            design, target, mol, top, ir_ = line.rstrip('\n').split('\t')
        except:
            print("CRITICAL: failed to parse line:\n%s\n" % line)
            sys.exit(1)
        ir = ir_.split(',')
        top = int(top)
        if design not in sites:
            sites[design] = {}
        if target not in sites[design]:
            sites[design][target] = {}
        if top not in sites[design][target]:
            sites[design][target][top] = {}
        sites[design][target][top][mol] = ir
    return sites
    
def read_contacts_map_from_file(file):
    contacts = {}
    if not os.path.isfile(file): return contacts
    lines = []
    with open(file, "r") as in_:
        fcntl.flock(in_, fcntl.LOCK_EX)
        lines = in_.readlines()
        fcntl.flock(in_, fcntl.LOCK_UN)
    for line in lines:
        try:
            design, target, top, conts_ = line.rstrip('\n').split('\t')
        except:
            print("CRITICAL: failed to parse line:\n%s\n" % line)
            sys.exit(1)
        conts = conts_.split(',')
        top = int(top)
        if design not in contacts:
            contacts[design] = {}
        if target not in contacts[design]:
            contacts[design][target] = {}
        contacts[design][target][top] = conts
    return contacts
    
def read_scores_from_file(file):
    scores = {}
    if not os.path.isfile(file): return scores
    lines = open(file).readlines()
    for line in lines:
        tmp = line.rstrip().split('\t')
        design, target, top, score, value = line.rstrip().split('\t')
        top = int(top)
        value = float(value)
        if design not in scores:
            scores[design] = {}
        if target not in scores[design]:
            scores[design][target] = {}
        if top not in scores[design][target]:
            scores[design][target][top] = {}
        scores[design][target][top][score] = value
    return scores
    
def read_masks_from_file(file):
    lines = open(file).readlines()
    masks = {}
    for line in lines:
        tmp = line.rstrip().split()
        try:
            masks[tmp[1]] = [ c for c in tmp[0] ]
        except:
            continue
        if 'seq' not in tmp[1]:
            masks[tmp[1]] = [ int(i) for i in masks[tmp[1]] ]
    return masks

def get_prediction_id(molid1,molid2,top):
    if top < 0:
        return MAX_PREDICTIONS + top
    return MAX_MOLECULES * MAX_MOLECULES * top + \
        MAX_MOLECULES * molid1 + \
        molid2
        
def assign_molecule_ids_from_file(file):
    lines = open(file).readlines()
    id = 0
    molids = {}
    for line in lines:
        tmp = line.split('\t')
        mol = tmp[0]
        if mol not in molids:
            molids[tmp[0]] = id
            id += 1
    return molids

    
build_folders()

DEFAULT_PARENT_UID = hash(DEFAULT_SEQUENCE)
DEFAULT_PARENT_SHORTUID = DEFAULT_PARENT_UID[0:N_SHORTUID_CHARS]
