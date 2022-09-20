from ase import Atoms, Atom 
from ase.build import bulk
from ase.visualize import view
from math import *
import numpy as np

# make components
l_cc = 1.35 # angstroms
l_cn = 1.47 # angstroms
l_ch = 1.09 # angstroms
l_co = 1.43
c_bond_dists = {
    "C":l_cc,
    "N":l_cn,
    "H":l_ch,
    "O":l_co
}

ext_a = pi/3
def make_carbon_ring(ring_atoms=["C","C","C","C","C","C"], connections=[]):
    
    c_ring = Atoms([Atom(ring_atoms[0], [0,0,0]),
                Atom(ring_atoms[1],[l_cc*sin(ext_a),-l_cc*cos(ext_a),0]),
                Atom(ring_atoms[2],[2*l_cc*sin(ext_a),0,0]),
                Atom(ring_atoms[3],[2*l_cc*sin(ext_a),l_cc,0]),
                Atom(ring_atoms[4],[l_cc*sin(ext_a),l_cc*(1+cos(ext_a)),0]),
                Atom(ring_atoms[5],[0,l_cc,0])])
    
    center = np.array([l_cc*sin(ext_a),l_cc/2,0])
    for i,element in enumerate(connections):
        if element != None:
            c_pos = c_ring[i].position
            dist_from_center = l_cc + c_bond_dists[element]
            rel_pos = np.subtract(c_pos, center)
            atom_dir = (rel_pos) / np.linalg.norm(rel_pos, axis=0)
            atom = Atom(element, (center+rel_pos)+(atom_dir*c_bond_dists[element]))
            c_ring.append(atom)
    return c_ring


def move_hydrogens_z(ring):
    ring = ring.copy()
    dz = .5
    for i,atom in enumerate(ring):
        if atom.symbol == 'H' and atom.position[1] == 1.895:
            atom.position[2] += dz
            dz *= -1
    return ring

def combine_two_rings_b(b1, b2):
    
    b1_c, b2_c = b1.copy(), b2.copy()
    b2_c.rotate(180, 'z')
    b2_c.translate([4*l_cc*sin(ext_a),0,0])
    b2_c.translate([l_cc*sin(pi/4),-l_cc*cos(pi/4),0])
    b1b2 = Atoms(cell=[4*l_cc*sin(ext_a)+2*l_cc*sin(pi/4),0,0])
    for atom in b1_c:
        b1b2.append(atom)
    for atom in b2_c:
        b1b2.append(atom)
    b1b2.center(vacuum=10, axis=(1,2))
    return b1b2


def combine_ring_simple_b(ring,simple):
    
    ring_c, simple_c = ring.copy(), simple.copy()
    ring_c = move_hydrogens_z(ring_c)    
    simple_c.rotate(180, 'z')
    simple_c.translate([2*l_cc*sin(ext_a),0,0])
    
    bond_len = c_bond_dists[simple_c[0].symbol]
    simple_c.translate([bond_len*sin(pi/4),-bond_len*cos(pi/4),0])
    
    b1b2 = Atoms(cell=[2*l_cc*sin(ext_a)+2*l_cc*sin(pi/4),0,0])
    for atom in ring_c:
        b1b2.append(atom)
    for atom in simple_c:
        b1b2.append(atom)
    b1b2.center(vacuum=10, axis=(1,2))
    return b1b2
    
#def combine_two_simple_b(b1, b2):
#    if b1[0].symbol == "C":
#        l_d = c_bond_dists[b2[0].symbol]
#    elif b1[0].symbol == "N":
#        if b2[0].symbol == "N":
#            l_d = 1.09
#        elif b2[0].symbol == "O":
#            l_d = 1.2
#    elif b1[0].symbol == "O":
#        if b2[0].symbol == "O":
#            l_d = 1.208
#    
#    b1_c, b2_c = b1.copy(), b2.copy()
#    b2_c.rotate(180, 'z')    
#    b2_c.translate([l_d*sin(pi/4),-l_d*cos(pi/4),0])
#    b1b2 = Atoms(cell=[l_d*sqrt(2),0,0])
#    for atom in b1_c:
#        b1b2.append(atom)
#    for atom in b2_c:
#        b1b2.append(atom)
#    b1b2.center(vacuum=10, axis=(1,2))
#    return b1b2

def combine_two_simple_b(b1, b2):
    if b1[0].symbol == "C":
        l_d = c_bond_dists[b2[0].symbol] + .15
    elif b1[0].symbol == "N":
        if b2[0].symbol == "N":
            l_d = 1.40
        elif b2[0].symbol == "O":
            l_d = 1.4
    elif b1[0].symbol == "O":
        if b2[0].symbol == "O":
            l_d = 1.408
    
    b1_c, b2_c = b1.copy(), b2.copy()
    b2_c.rotate(180, 'z')    
    b2_c.translate([l_d*sin(1.5*pi/4),-l_d*cos(1.5*pi/4),0])
    b1b2 = Atoms(cell=[2*l_d*sin(1.5*pi/4),0,0])
    for atom in b1_c:
        b1b2.append(atom)
    for atom in b2_c:
        b1b2.append(atom)

    b1b2.center(vacuum=15, axis=(1,2))
    return b1b2
     

def combine_blocks_b(block1, block2, length=2):
    
    if len(block1) >= 6 and len(block2) >= 6: # if we have two rings
        return combine_two_rings_b(block1, block2)
    elif len(block1) >= 6 or len(block2) >= 6: ### for now just the 1 or 2-atom blocks
        return combine_ring_simple_b(block1,block2)
    else:
        return combine_two_simple_b(block1,block2)




########################################################
########## The rest is for 'C' type molecules ##########
############## with a hexagonal structure ##############
########################################################
def combine_two_rings_c(b1, b2):
    if 'B' in b1.symbols and 'B' in b2.symbols:
        l_d = 1.94 # distance between blocks!
    elif 'B' in b1.symbols or 'B' in b2.symbols:
        l_d = 1.56
    else:
        l_d = 1.35
    b1_c = b1.copy()
    b2_c1,b2_c2,b2_c3 = b2.copy(), b2.copy(), b2.copy()
    
#    b2_c1.rotate(60,'z',center=[l_cc*sin(ext_a),l_cc/2,0])
#    b2_c1.translate([2*l_cc*sin(ext_a),-l_cc,0])
#    b2_c1.translate([l_d*sin(pi/3),-l_d*cos(pi/3),0])
    
#    b2_c2.rotate(60,'z',center=[l_cc*sin(ext_a),l_cc/2,0])
#    b2_c2.translate([-2*l_cc*sin(ext_a),-l_cc,0])
#    b2_c2.translate([-l_d*sin(pi/3),-l_d*cos(pi/3),0])

    b2_c3.translate([0,2*l_cc+l_d,0])
    b2_c3.rotate(60,'z',center=[l_cc*sin(ext_a),(2*l_cc+l_d)+l_cc/2,0]) 
    
    
    xvec = [6*l_cc*sin(pi/3),0,0]
    yvec = [2*l_cc*sin(pi/3)+(l_d*sin(pi/4)),
            4*l_cc+(l_d*sin(pi/4)),
            0]
    zvec = [0,0,0]
    b1b2 = Atoms(cell=[xvec,yvec,zvec])
    for atom in b1_c:
        b1b2.append(atom)
#     for atom in b2_c1:
#         b1b2.append(atom)
#     for atom in b2_c2:
#         b1b2.append(atom)
    for atom in b2_c3:
        b1b2.append(atom)
        
    b1b2.center(vacuum=10, axis=(2))
    return b1b2


def combine_ring_simple_c(ring,simple):
    if 'B' in ring.symbols:
        if simple[0].symbol == 'N':
            l_d = 1.752
        elif simple[0].symbol == 'C':
            l_d = 1.56
    else:
        l_d = c_bond_dists[simple[0].symbol]
    ring_c = ring.copy()
    simple_c1,simple_c2,simple_c3 = simple.copy(), simple.copy(), simple.copy()
    
    simple_c1.translate([2*l_cc*sin(ext_a),0,0])
    simple_c1.translate([l_d*sin(pi/3),-l_d*cos(pi/3),0])
    
    simple_c2.translate([-l_d*sin(pi/3),-l_d*cos(pi/3),0])

    simple_c3.translate([0,2*l_cc+l_d,0])
    simple_c3.rotate(60,'z',center=[l_cc*sin(ext_a),(2*l_cc+l_d)+l_cc/2,0]) 
    b1b2 = Atoms()
    for atom in ring_c:
        b1b2.append(atom)
    for atom in simple_c1:
        b1b2.append(atom)
    for atom in simple_c2:
        b1b2.append(atom)
    for atom in simple_c3:
        b1b2.append(atom)
    return b1b2 
    
def combine_two_simple_c(b1, b2):
    if b1[0].symbol == "N" and b2[0].symbol == "N":
        l_d = 1.3
    elif b1[0].symbol == "C":
        l_d = c_bond_dists[b2[0].symbol]
    elif b2[0].symbol == "C":
        l_d = c_bond_dists[b1[0].symbol]
        
    b1_c = b1.copy()
    b2_c1,b2_c2,b2_c3 = b2.copy(), b2.copy(), b2.copy()    
    
    b2_c1.translate([l_d*sin(pi/3),-l_d*cos(pi/3),0])
    
    b2_c2.translate([-l_d*sin(pi/3),-l_d*cos(pi/3),0])
    
    b2_c3.translate([0,l_d,0])
    
    
    xvec = [2*l_d*sin(pi/3),0,0]
    yvec = [l_d*sin(pi/3),(3/2)*l_d,0]
    zvec = [0,0,0]
    b1b2 = Atoms(cell=[xvec,yvec,zvec])
    for atom in b1_c:
        b1b2.append(atom)
#     for atom in b2_c1:
#         b1b2.append(atom)
#     for atom in b2_c2:
#         b1b2.append(atom)
    for atom in b2_c3:
        b1b2.append(atom)
        
    b1b2.center(vacuum=10, axis=2)
    return b1b2
    
    
def combine_blocks_c(block1, block2):
    if len(block1) >= 6 and len(block2) >= 6: # if we have two rings
        return combine_two_rings_c(block1, block2)
    elif len(block1) >= 6 or len(block2) >= 6: ### for now just the 1 or 2-atom blocks
        return combine_ring_simple_c(block1,block2)
    else:
        return combine_two_simple_c(block1,block2)
