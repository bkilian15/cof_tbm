def set_control_parameters(calculation=None, restart_mode=None,
                          prefix=None, tstress=None, tprnfor=None, etot_conv_thr=None,
                           forc_conv_thr=None,pseudo_dir=None, outdir=None, verbosity=None):
    control = {}    

    control["calculation"] = calculation
    if restart_mode != None:
        control["restart_mode"] = restart_mode
    control["prefix"] = prefix
    if tstress != None:
        control["tstress"] = tstress
    if tprnfor != None:
        control["tprnfor"] = tprnfor
    if etot_conv_thr != None:
        control["etot_conv_thr"] = etot_conv_thr
    if forc_conv_thr != None:
        control["forc_conv_thr"] = forc_conv_thr
    if pseudo_dir != None:
        control["pseudo_dir"] = pseudo_dir
    else:
        control["pseudo_dir"] = f"./{prefix}/pseudo_dir"
    if outdir != None:
        control["outdir"] = outdir
    else:
        control["outdir"] = f"./{prefix}/outdir"
    if verbosity != None:
        control["verbosity"] = verbosity
    
    return control

def set_system_parameters(ibrav=None, celldm=None, nat=None, ntyp=None,
                            nspin=None,  starting_magnetization=None,
                            ecutwfc=None, ecutrho=None, occupations=None,
                            smearing=None, degauss=None, nbnd=None):
    system = {}
    if ibrav != None: system["ibrav"] = ibrav
    if celldm != None: system["celldm"] = celldm
    if nat != None: system["nat"] = nat
    if ntyp != None: system["ntyp"] = ntyp
    if nspin != None: system["nspin"] = nspin
    if starting_magnetization != None: 
        system["starting_magnetization"] = starting_magnetization
    if ecutwfc != None: system["ecutwfc"] = ecutwfc
    if ecutrho != None: system["ecutrho"] = ecutrho
    if occupations != None: system["occupations"] = occupations
    if smearing != None: system["smearing"] = smearing
    if degauss != None: system["degauss"] = degauss
    if nbnd != None: system["nbnd"] = nbnd

    return system

def set_electrons_parameters(electron_maxstep=None,diagonalization=None,
                            mixing_mode=None,mixing_beta=None,
                            conv_thr=None):
    electrons = {}
    if electron_maxstep != None: electrons["electron_maxstep"] = electron_maxstep
    if diagonalization != None: electrons["diagonalization"] = diagonalization
    if mixing_mode != None: electrons["mixing_mode"] = mixing_mode
    if mixing_beta != None: electrons["mixing_beta"] = mixing_beta
    if conv_thr != None: electrons["conv_thr"] = conv_thr
        
    return electrons


def set_cell_parameters(cell_dynamics=None, press=None, wmass=None,
                       cell_factor=None,press_conv_thr=None,cell_dofree=None):

    cell = {}

    if cell_dynamics != None: cell["cell_dynamics"] = cell_dynamics
    if press != None: cell["press"] = press
    if wmass != None: cell["wmass"] = wmass
    if cell_factor != None: cell["cell_factor"] = cell_factor
    if press_conv_thr != None: cell["press_conv_thr"] = press_conv_thr
    if cell_dofree != None: cell["cell_dofree"] = cell_dofree
   
    return cell
