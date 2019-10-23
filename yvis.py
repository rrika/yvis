import os.path
import array
import argparse
import subprocess
import sys
sys.path.append("srctools")
import srctools.bsp
from srctools.bsp import BSP_LUMPS

def nostalgia_parser():
	parser = argparse.ArgumentParser(usage='%(prog)s [-fast] [-game DIRECTORY] bspfile', add_help=False)
	common_opts = parser.add_argument_group('Common options')
	common_opts.add_argument('bspfile')
	common_opts.add_argument('-h', '--help', action='help',
		help=argparse.SUPPRESS)
	common_opts.add_argument('-v', '-verbose', dest='verbose', action='store_true',
		help='Turn on verbose output.')
	common_opts.add_argument('-fast', dest='fast', action='store_true',
		help='Only do first quick pass on vis calculations.')
	common_opts.add_argument('-low', dest='low', action='store_true',
		help='Run as an idle-priority process.')
	common_opts.add_argument('-game', '-vproject', metavar='<directory>',
		help='Override the VPROJECT environment variable.')

	other_opts = parser.add_argument_group('Other options')
	other_opts.add_argument('-radius_override', metavar='<n>', dest='radius',
		help='Force a vis radius, regardless of whether an env_fog_controller specifies one.')
	other_opts.add_argument('-threads', metavar='<n>', dest='nthreads',
		help='Control the number of threads vbsp uses (defaults to the # of processors on your machine).')
	other_opts.add_argument('-trace', metavar=('<start cluster>', '<end cluster>'), nargs=2, dest='trace',
		help='Writes a linefile that traces the vis from one cluster to another for debugging map vis.')

	other_opts.add_argument('-FullMinidumps',        action='store_true', help=argparse.SUPPRESS)
	other_opts.add_argument('-nox360', '-nosort',    action='store_true', help=argparse.SUPPRESS)
	other_opts.add_argument('-mpi', '-novconfig',    action='store_true', help=argparse.SUPPRESS)
	other_opts.add_argument('-tmpin', '-tmpout',     action='store_true', help=argparse.SUPPRESS)
	other_opts.add_argument('-allowdebug', '-steam', action='store_true', help=argparse.SUPPRESS)
	other_opts.add_argument('-mpi_pw',               nargs=1, help=argparse.SUPPRESS)
	return parser

def plain_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('bspfile')
	parser.add_argument('--prt', metavar='<prtfile>', dest='prtfile')
	parser.add_argument('--game', metavar='<directory>', dest='game')
	return parser


def determine_radius(vmf):
	for entity in vmf.entities:
		if entity["classname"].casefold() == "env_fog_controller":
			try:
				return float(entity.get("farz", ""))
			except ValueError:
				continue

def build_cluster_table():
	clusters = [[] for i in range(numclusters)]
	for i, leaf in enumerate(dleafs):
		clusters[leaf.cluster].append(i)
	return clusters

def iter_leaf_leaf_vis():
	pass

#def iter_cluster_faces():
#	firstleafface = dleafs[nClusterLeaf].firstleafface;
#	numleaffaces = dleafs[nClusterLeaf].numleaffaces
#	for leafFaceID in range(firstleafface, firstleafface + numleaffaces):
#		int faceID = dleaffaces[nFirstFaceID + leafFaceID];

# contents, cluster, area_flags, mins, maxs, firstleafface, numleaffaces, firstleafbrush, numleafbrushes, leafWaterDataID
dleaf_1_fmt = "<ihh3h3hHHHHh"

CONTENTS_SOLID         =  0x001
CONTENTS_SLIME         =  0x010
CONTENTS_TESTFOGVOLUME =  0x100
SURF_WARP              = 0x0008

def calc_visible_fog_volumes(dleaf, iter_leaf_leaf_vis):
	# for every water leaf i:
	#    for every air leaf j:
	#       dleafs[j].contents |= CONTENTS_TESTFOGVOLUME

	for i in range(numleafs):
		contents = dleafs[i].contents
		waterid  = dleafs[i].leafWaterDataID

		if contents & (CONTENTS_SOLID | CONTENTS_SLIME):
			continue
		if waterid == -1:
			continue

		for j in iter_leaf_leaf_vis(i, DVIS_PVS):
			jcontents = dleafs[j].contents
			jwaterid  = dleafs[j].leafWaterDataID

			if jcontents & CONTENTS_SOLID:
				continue
			if jwaterid != -1:
				continue

			dleafs[nClusterLeaf].contents |= CONTENTS_TESTFOGVOLUME
"""
def CalcDistanceFromLeafToWater(leafNum):
	int j, k;

	# If we know that this one doesn't see a water surface then don't bother doing anything.
	if ((dleafs[leafNum].contents & CONTENTS_TESTFOGVOLUME) == 0) && ( dleafs[leafNum].leafWaterDataID == -1 ):
		return 65535
	
	# First get the vis data..
	int cluster = dleafs[leafNum].cluster;
	if (cluster < 0)
		return 65535; # FIXME: make a define for this.
	
	float minDist = 65535.0f; # FIXME: make a define for this.
	
	Vector leafMin, leafMax;
	
	leafMin[0] = ( float )dleafs[leafNum].mins[0];
	leafMin[1] = ( float )dleafs[leafNum].mins[1];
	leafMin[2] = ( float )dleafs[leafNum].mins[2];
	leafMax[0] = ( float )dleafs[leafNum].maxs[0];
	leafMax[1] = ( float )dleafs[leafNum].maxs[1];
	leafMax[2] = ( float )dleafs[leafNum].maxs[2];

	# something something build convex hull
	for j in iter_leaf_leaf_vis(leafNum, DVIS_PVS):
		jcontents = dleaf[j].contents
		jwaterid  = dleaf[j].leafWaterDataID

		if (jcontents & CONTENTS_TESTFOGVOLUME) == 0 && jwaterid == -1:
			continue

		for faceId in iter_cluster_faces(j):
			dface_t *pFace = &dfaces[faceID]
			if( pFace->texinfo == -1 )
				continue;

			texinfo_t *pTexInfo = &texinfo[pFace->texinfo];
			if( pTexInfo->flags & SURF_WARP )
				# Woo hoo!!!  We found a water face.
				# compare the bounding box of the face with the bounding
				# box of the leaf that we are looking from and see
				# what the closest distance is.
				# FIXME: this could be a face/face distance between the water
				# face and the bounding volume of the leaf.
				
				# Get the bounding box of the face
				Vector faceMin, faceMax;
				GetBoundsForFace( faceID, faceMin, faceMax );
				float dist = GetMinDistanceBetweenBoundingBoxes( leafMin, leafMax, faceMin, faceMax );
				if( dist < minDist )
					minDist = dist;

	return minDist
"""

def calc_distance_from_leaves_to_water():
	a = array.array('H', [calc_distance_from_leaf_to_water(i) for i in range(numleafs)])
	if sys.byteorder != 'little':
		a.byteswap()
	return a


def main():
	args = nostalgia_parser().parse_args()
	bsp = srctools.bsp.BSP(args.bspfile, srctools.bsp.VERSIONS.PORTAL)

	if not hasattr(args, "prtfile") or not args.prtfile:
		if ".bsp" in args.bspfile:
			args.prtfile = args.bspfile.replace(".bsp", ".prt")
		else:
			args.prtfile = args.bspfile + ".prt"

	# with open(prtfile, "r") as f:
	# 	prt = f.read()
	# nclusters = int(prt.split(" ", 3)[1])

	# to read cluster assignment
	lump_leafs     = bsp.get_lump(BSP_LUMPS.LEAFS)

	# to measure distance to water faces
	lump_faces     = bsp.get_lump(BSP_LUMPS.FACES)
	lump_leaffaces = bsp.get_lump(BSP_LUMPS.LEAFFACES)
	lump_edges     = bsp.get_lump(BSP_LUMPS.EDGES)
	lump_surfedges = bsp.get_lump(BSP_LUMPS.SURFEDGES)
	lump_vertexes  = bsp.get_lump(BSP_LUMPS.VERTEXES)
	lump_texinfo   = bsp.get_lump(BSP_LUMPS.TEXINFO)

	# to extract a vis radius from an env_fog_controller
	vmf = bsp.read_ent_data()

	subprocess.check_call(["cargo", "run", "--release", "--manifest-path=yvis/Cargo.toml", args.prtfile])

	build_cluster_table()

	lump_leafs     = calc_visible_fog_volumes(lump_leafs)
	lump_waterdist = calc_distance_from_leaves_to_water()

	bsp.get_lump(BSP_LUMPS.LEAFS).data = lump_leafs
	bsp.get_lump(BSP_LUMPS.VISIBILITY).data = lump_visibility
	bsp.get_lump(BSP_LUMPS.LEAFMINDISTTOWATER).data = lump_waterdist

	bsp.save(os.path.basename(args.bspfile)+".withvis")

if __name__ == '__main__':
	main()
