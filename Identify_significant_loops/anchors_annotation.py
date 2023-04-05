##### 
##### HiChIP Loops anchors annotation
##### author: Zhen Y

##### usage: anchors_annotation.py -1 left_anchor_enhancer.bed -2 right_anchor_enhancer.bed -3 left_anchor_promoter.bed -4
##### left_anchor_promoter.bed -o output.bed
##### anchors of each HiChIP loop were annotate with enhancers and promoters with bedtools intersect function
##### left_anchor_enhancer.bed: left anchors overlap with enhancers
##### right_anchor_enhancer.bed: right anchors overlap with enhancers
##### left_anchor_promoter.bed: left anchors overlap with promoters
##### right_anchor_promoter.bed: right anchors overlap with promoters
##### output.bed: output file name with loops annotation (EE, EP, PE, PP)

import argparse


def main_annotation(l_e, r_e, l_p, r_p, output):

    left_coords = {}
    right_coords_p = {}
    right_coords_e = {}
    loops = {}
    
    with open(l_e, 'r') as f:
        for line in f:
            entry = line.strip().split() # read each line from left_enhancer anchors
            id = entry[3]  # file format: chr start end loop_id (tab delimited)
            # create loops dictionary with loop_id as key            
            loops[id] = {}
            loops[id]['left'] = ['E']
            
            left_coords[id] = [entry[0], entry[1], entry[2], entry[4]] # coordinates + enhancer ID
               
    with open(r_e, 'r') as f:
        for line in f:
            entry = line.strip().split()  # read each line from right_enhancer anchors
            id = entry[3]
            
            if id in loops:
                loops[id]['right'] = ['E']
            else:
                loops[id] = {}
                loops[id]['right'] = ['E']
    
            right_coords_e[id] = [entry[0], entry[1], entry[2], entry[4], entry[5]] # coordinates + enhancer_id + reads support
            
    with open(l_p, 'r') as f:
        for line in f:
            entry = line.strip().split()  # read each line from left_promoter anchors
            id = entry[3]
            if id in loops:
                if 'left' in loops[id]:
                    loops[id]['left'].append('P')
                else:
                    loops[id]['left'] = ['P']
            else:
                loops[id] = {}
                loops[id]['left'] = ['P']
            
            left_coords[id] = [entry[0], entry[1], entry[2], entry[4]]  # coordinates + genes

    with open(r_p, 'r') as f:
        for line in f:
            entry = line.strip().split()  # read each line from right_promoter anchors
            id = entry[3]
            
            if id in loops:
                if 'right' in loops[id]:
                    loops[id]['right'].append('P')
                else:
                    loops[id]['right'] = ['P']
            else:
                loops[id] = {}
                loops[id]['right'] = ['P']
            
            right_coords_p[id] = [entry[0], entry[1], entry[2], entry[4], entry[5]] # coordinates + genes + reads support


    with open(output, 'w') as out:   # saving loops and annotation to file
        for id in sorted(loops.keys()):
            if 'left' not in loops[id] or 'right' not in loops[id]:
                continue
            
            if 'P' in loops[id]['left'] and 'P' in loops[id]['right']:
                annot = 'PP'
                tdf = left_coords[id] + right_coords_p[id] + [annot]
                out.write('\t'.join(tdf) + '\n')

            if 'E' in loops[id]['left'] and 'P' in loops[id]['right']:
                annot = 'EP'
                tdf = left_coords[id] + right_coords_p[id] + [annot]
                out.write('\t'.join(tdf) + '\n')

            if 'P' in loops[id]['left'] and 'E' in loops[id]['right']:
                annot = 'PE'
                tdf = left_coords[id] + right_coords_e[id] + [annot]
                out.write('\t'.join(tdf) + '\n')

            if 'E' in loops[id]['left'] and 'E' in loops[id]['right']:
                annot = 'EE'
                tdf = left_coords[id] + right_coords_e[id] + [annot]
                out.write('\t'.join(tdf) + '\n')

            



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Loops annotation')
    parser.add_argument('-1', '--left_enhancer', metavar='left enhancer', required=True, 
        help='input left_anchor_enhancer bed file')
    parser.add_argument('-2', '--right_enhancer', metavar='right enhancer', required=True, 
        help='input right_anchor_enhancer bed file')
    parser.add_argument('-3', '--left_promoter', metavar='left promoter', required=True, 
        help='input left_anchor_promoter bed file')
    parser.add_argument('-4', '--right_promoter', metavar='right promoter', required=True, 
        help='input right_anchor_promoter bed file')
    parser.add_argument('-o', '--output', metavar='output file', required=True, 
        help='direct annotation loops to file')
    args = parser.parse_args()

    l_e = args.left_enhancer
    r_e = args.right_enhancer
    l_p = args.left_promoter
    r_p = args.right_promoter
    output = args.output

    main_annotation(l_e, r_e, l_p, r_p, output)