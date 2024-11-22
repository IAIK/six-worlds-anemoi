from os import getcwd, walk, path, makedirs
import re
import sys

load('complexities.sage')

def parseFileName(filePath):
    pattern = r'anemoi_p(?P<p>\d+|[A-Z]+\d*-\d+)_n(?P<n>\d+)_l(?P<l>\d+)_a(?P<a>\d+)_m(?P<m>[A-Z]+)(?P<noll>_noll)?(?:_o(?P<o>\d+))?'
    match = re.search(pattern, filePath)

    if match:
        p_ = match.group("p")
        p = int(p_) if (len(p_.split("-")) == 1) else p_
        n = int(match.group("n"))
        l = int(match.group("l"))
        a = int(match.group("a"))
        m = match.group("m")
        noll = False if match.group("noll") else True
        o = int(match.group("o")) if match.group("o") else 1  # Default value 1 if no ordering is specified
    else:
        print(f"No match: {filePath}")

    return p, n, l, a, m, noll, o


def parseOutputFile(filePath):
    
    try:
        with open(filePath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError as e:
        print(e)
        return []
    
    attackInfo = []
    print(filePath)
    
    for i, line in enumerate(lines):
        if line.startswith("n_r:"):
            pattern = r"n_r: (?P<n_r>\d+), n_v: (?P<n_v>\d+), n_e: (?P<n_e>\d+), mdeg: (?P<mdeg>\d+)"
            match = re.search(pattern, line)

            #if match:
            n_r = match.group("n_r")
            n_v = match.group("n_v")
            n_e = match.group("n_e")
            mdeg = match.group("mdeg")
            
            attackInfo.append({'nr': Integer(n_r), 'nv': Integer(n_v), 'ne': Integer(n_e), 'mdeg': Integer(mdeg), 'degs': [None],
                               'GB': {'t': -1, 'dregs': [None], 'dreg': None, 'c': {2: None, 2.37: None, 3: None}}, # Step 1: GB
                               'FGLM': {'t': -1, 'vsdim': None, 'c': {2: None, 2.37: None, 3: None}},               # Step 2: FGLM
                               'ELIM': {'t': -1, 'unideg': None, 'nsols': None, 'nsols_mult': None, 'c': None},     # Step 3: ELIM / Variety
                               'crashed': False
                              })   
        
        elif line.startswith("Magma crashed"):
            attackInfo[-1]['crashed'] = True
            
        elif line.startswith("degs = "): # Degrees of polynomials
            degs = [Integer(x) for x in line[len("degs = ")+1:-2].split(', ')]
            attackInfo[-1]['degs'] = degs
            
        elif line.startswith("Get rep mat (dim: "):
            vsdim = Integer(line[len("Get rep mat (dim: "):-2])
            attackInfo[-1]['FGLM']['vsdim'] = vsdim
        
        elif line.startswith("Quotient dimension: "):
            vsdim = Integer(line[len("Quotient dimension: "):-1])
            attackInfo[-1]['FGLM']['vsdim'] = vsdim
        
        elif line.startswith("nsols_mult ="): # ELIM / Variety
            nsols = Integer(line[len("nsols_mult ="):-1])
            attackInfo[-1]['ELIM']['nsols_mult'] = nsols
        
        #elif line.startswith("Timing (ELIM): ") and attackInfo[-1]['t3'] == -1:
        #    t = float(line[len("Timing (ELIM): "):-2])
        #    attackInfo[-1]['t3'] = t
        
        elif line.startswith("Time: "):
            t = line[len("Time: "):-1]
            t = float(t.strip('[r]'))
            
            next_line = lines[i+1]
            if next_line.startswith("dregs = "): # GB
                # dregs can be empty if immediately done
                dregs = [Integer(x) for x in next_line[len("dregs = ")+1:-2].split(', ') if x != '']
                if len(dregs) == 0:
                    print("--> GB finished without performing any steps.")
                else:
                    attackInfo[-1]['GB']['dregs'] = dregs
                    attackInfo[-1]['GB']['dreg'] = max(dregs)
                    attackInfo[-1]['GB']['t'] = t
                    for w in OMEGAS:
                         attackInfo[-1]['GB']['c'][w] = gb_complexity(ne=attackInfo[-1]['ne'], nv=attackInfo[-1]['nv'], 
                                                                      dreg=attackInfo[-1]['GB']['dreg'], degs=attackInfo[-1]['degs'], w=w)
            elif next_line.startswith("unideg ="): # FGLM
                unideg = Integer(next_line[len("unideg = "):-1])
                
                if attackInfo[-1]['FGLM']['vsdim'] == None:
                    # This happens if after GB computation, GB is already lex GB  -> for F2n, a=3
                    attackInfo[-1]['FGLM']['vsdim'] = unideg
                    print("--> FGLM finished without performing any steps.")
                assert(attackInfo[-1]['FGLM']['vsdim'] == unideg) # if this is the case, then shape position
                
                attackInfo[-1]['ELIM']['unideg'] = unideg
                attackInfo[-1]['FGLM']['t'] = t
                for w in OMEGAS:
                    attackInfo[-1]['FGLM']['c'][w] = fglm_complexity(nv=attackInfo[-1]['nv'], dI=attackInfo[-1]['FGLM']['vsdim'], w=w)
            elif next_line.startswith("nsols ="): # ELIM / Variety
                nsols = Integer(next_line[len("nsols ="):-1])
                attackInfo[-1]['ELIM']['nsols'] = nsols
                attackInfo[-1]['ELIM']['t'] = t
                c_fac = fac_complexity(d=attackInfo[-1]['ELIM']['unideg'])
                attackInfo[-1]['ELIM']['c'] = c_fac
            
    return attackInfo


if __name__ == "__main__":
    args = sys.argv
    print(args)
    if len(args) == 2: # parse single file
        filePath = args[1]
        results = parseOutputFile(filePath)
        save(results,filePath.replace("/","_"))
        
    else: # parse all files in result directory
        resultFiles = []
        dirname = "results/"
        for (dirPath, dirNames, fileNames) in walk(dirname):
            resultFiles += [dirPath + fileName for fileName in fileNames]
            break
    
        # Prime field results: {p: {a : {m : {o : attackInfo}}}}
        # Binary field results: {n: {a : {m : {o : attackInfo}}}}
        results = {}
        results_Fp = {}
        results_F2n = {}

        for filePath in resultFiles:
            p, n, l, a, m, noll, o = parseFileName(filePath)
            results = results_Fp if n == 1 else results_F2n
            key = (p,n)

            if key not in results:
                results[key] = {}
            if a not in results[key]:
                results[key][a] = {}
            if m not in results[key][a]:
                results[key][a][m] = {}

            results[key][a][m][o] = parseOutputFile(filePath)

        save(results_Fp,'results_Fp')
        save(results_F2n,'results_F2n')