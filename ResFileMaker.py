#! /usr/bin/env python


import optparse    # for option sorting
import rosetta.core.pack.task    # for using resfiles

def generate_resfile_from_pose(pose, resfilename,
        pack = True, design = False, input_sc = True,
        freeze = [], specific = {}):
    """
    Writes a resfile for  <pose>  named  <resfilename>
       <pack> = True allows packing by default
       <design> = True allows design using all amino acids by default
       <input_sc> = True allows usage of the original side chain conformation
       <freeze> is an optional list of (pose) residue numbers to exclude
            (preserve the side chain conformations of these residues)
       <specific> is an optional dictionary with (pose) residue numbers as keys
            and resfile keywords as corresponding values
            (for setting individual residue options, it may be easier to add
            these numbers to freeze and edit the resfile manually)

    example:
        generate_resfile_from_pose(pose,'1YY8.resfile')
    See also:
        Pose
        PackRotamersMover
        TaskFactory

    """
    # determine the header, default settings
    header = ''
    if pack:
        if not design:
            header += 'NATAA\n'
        else:
            header += 'ALLAA\n# ALLAA will NOT work on bridged Cysteines\n'
    else:
        header += 'NATRO\n'
    if input_sc:
        header += 'USE_INPUT_SC\n'
    to_write = header + 'start\n'
    # add  <freeze>  list to  <specific>  dict
    for i in freeze:
        specific[i] = 'NATRO'
    #  <specific>  is a dictionary with keys() as pose resi numbers
    #    and values as resfile keywords (PIKAA
    # use PDBInfo object to write the resfile
    info = pose.pdb_info()
    # pose_from_sequence returns empty PDBInfo, Pose() makes NULL
    if info and info.nres():
        for i in specific.keys():
            num = pose.pdb_info().number(i)
            chain = pose.pdb_info().chain(i)
            to_write += str(num).rjust(4) + str(chain).rjust(3) + '  ' + specific[i] + '  \n'
    else:
        for i in specific.keys():
            num = i
            chain = ' '
            to_write += str(num).rjust(4) + str(chain).rjust(3) + '  ' + specific[i] + '  \n'
    f = open(resfilename,'w')
    f.write(to_write)
    f.close()