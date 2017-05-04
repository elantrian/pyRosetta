# !/usr/bin/python
#:noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available
# (c) under license.
# (c) The Rosetta software is developed by the contributing members of the
# (c) Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington UW
# (c) TechTransfer, email: license@u.washington.edu.
"""Brief:   This PyRosetta script mutates residues and calculates delta rosetta score


Author: Elizabeth Lagesse, with excerpts borrowed from Jason Labonte's
delta_score_per_mutation.py

"""

import sys
import argparse
from rosetta import init, pose_from_pdb, get_fa_scorefxn, \
                    standard_packer_task, change_cys_state, \
                    Pose, MoveMap, RotamerTrialsMover, MinMover
from toolbox import mutate_residue

PACK_RADIUS = 10.0

mutationsList = ['7S', '19E']

#make a dictionary of mutations to make
mutationDict = {} #key: position, value: new AA to mutate too

for element in mutationsList:
	thisLetter = element[-1]
	mutationDict[element]=thisLetter

print mutationDict

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('pdb_filename',
                    help='the filename of the PDB structure to be evaluated')
parser.add_argument('-m', '--minimize', action='store_true',
                    help='flag to perform minimization after each mutation')
args = parser.parse_args()

init(extra_options='-mute basic -mute core')

data = ['Variant,Rosetta Score,"delta-delta-G"\n']

initial_pose = pose_from_pdb(args.pdb_filename)
sf = get_fa_scorefxn()
mm = MoveMap()
mm.set_bb(True)
mm.set_chi(True)
pre_pre_packing_score = sf(initial_pose)

task = standard_packer_task(initial_pose)
task.restrict_to_repacking()
task.or_include_current(True)
pack_rotamers_mover = RotamerTrialsMover(sf, task)
pack_rotamers_mover.apply(initial_pose)

min_mover = MinMover()
min_mover.movemap(mm)
min_mover.score_function(sf)
min_mover.min_type('linmin')
if args.minimize:
    min_mover.apply(initial_pose)
    initial_pose.dump_pdb(str("SUMO_WT_min.pdb"))

post_pre_packing_score = sf(initial_pose)

print
print 'Reference Protein:', args.pdb_filename
print '  Score:'
print '    Before pre-packing:', pre_pre_packing_score
print '    After pre-packing:', post_pre_packing_score
print

data.append('WT,' + str(post_pre_packing_score) + ',0.0\n')

print 'Making mutations',
if args.minimize:
	print 'with minimization',
print 'and scoring poses. ',

#mutate residue and minimize
for entry in sorted(mutationDict):
	seq_pos = int(entry[:-1])
	AA = mutationDict[entry]
	variant_name = str(entry)
	mutant_pose = mutate_residue(initial_pose, seq_pos, AA,
                                             PACK_RADIUS, sf)
	if args.minimize:
		min_mover.apply(mutant_pose)
		mutant_pose.dump_pdb(str("SUMO_"+variant_name+"_min.pdb"))
		variant_score = sf(mutant_pose)

		data.append(variant_name + "," + \
                str(variant_score) + "," + \
                str(variant_score - post_pre_packing_score) + "\n")


# Output results.
data_filename = args.pdb_filename[:-4] + '_variant_scores.csv'
with open(data_filename, "w") as f:
    f.writelines(data)

print 'Data written to:', data_filename

