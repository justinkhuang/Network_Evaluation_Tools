###################################################################
# Command line script to analyze network on node sets of interest #
###################################################################

from network_evaluation_tools import network_evaluation_functions as nef
from network_evaluation_tools import data_import_tools as dit
from network_evaluation_tools import gene_conversion_tools as gct
import argparse
import os
import pandas as pd

# Checking valid alpha and p values (Range is 0.0-1.0 exclusive)
# Value can also be None.
def restricted_float(x):
	if x is not None:
		x = float(x)
		if x <= 0.0 or x >= 1.0:
			raise argparse.ArgumentTypeError("%r not in range (0.0, 1.0) exclusive"%(x,))
	return x

# Checking valid integer values (for all values that must be >0)
def positive_int(x):
	x = int(x)
	if x <= 0:
		 raise argparse.ArgumentTypeError("%s must be a positive integer" % x)
	return x

# Valid file path check (Does not check file formatting, but checks if given path exists and is readable)
def valid_infile(in_file):
	if not os.path.isfile(in_file):
		raise argparse.ArgumentTypeError("{0} is not a valid input file path".format(in_file))	
	if os.access(in_file, os.R_OK):
		return in_file
	else:
		raise argparse.ArgumentTypeError("{0} is not a readable input file".format(in_file))

# Valid output directory path check (Checks if the output directory path can be found and written to by removing given filename from full path)
# Note: This uses '/' character for splitting pathnames on Linux and Mac OSX. The character may need to be changed to '\' for Windows executions
def valid_outfile(out_file):
	outdir = '/'.join(out_file.split('/')[:-1])
	if not os.path.isdir(outdir):
		raise argparse.ArgumentTypeError("{0} is not a valid output directory".format(outdir))
	if os.access(outdir, os.W_OK):
		return out_file
	else:
		raise argparse.ArgumentTypeError("{0} is not a writable output directory".format(outdir))

if __name__ == "__main__":
	# Network Evaluation Setup Variables
	parser = argparse.ArgumentParser(description='Analyze network performance on ability to aggregate sets of nodes in network space.')
	parser.add_argument("network_path", type=valid_infile, 
		help='Path to file of network to be evaluated. File must be 2-column edge list where each line is a gene interaction separated by a common delimiter.')
	parser.add_argument("node_sets_file", type=valid_infile, 
		help='Path to file of node sets. Each line is a list, separated by a common delimiter. The first item in each line will be the name of the node set.')
	parser.add_argument("actual_AUPRCs_save_path", type=valid_outfile, 
		help='CSV file path of network evaluation result scores (AUPRCs). This script minimally returns these values to save. Must have a writable directory.')		
	parser.add_argument('-v', '--verbose', default=False, action="store_true", required=False,
		help='Verbosity flag for reporting on patient similarity network construction steps.')	
	parser.add_argument('-netd', '--net_file_delim', type=str, default='\t', required=False,
		help='Delimiter used in network file between columns. Default is tab white space.')
	parser.add_argument('-setd', '--set_file_delim', type=str, default='\t', required=False,
		help='Delimiter used in node set file to delimit lists. Default is tab white space.')	
	parser.add_argument("-p", "--sample_p", type=restricted_float, default=None, required=False,
		help='Sub-sampling percentage for node sets of interest. Default is None. Each gene set''s p is automatically determined by the network in this case.')
	parser.add_argument("-a", "--alpha", type=restricted_float, default=None, required=False,
		help='Propagation constant to use in the propagation of node sub-samples over given network. Overrides alpha calculation model if given.')
	parser.add_argument("-n", "--sub_sample_iter", type=positive_int, default=30, required=False,
		help='Number of times to perform sub-sampling during performance recovery (AUPRC) calculation for each node set. Default is 30.')
	parser.add_argument('-c', '--cores', type=positive_int, default=1, required=False,
		help='Number of cores to be utilized by machine for performance calculation step. NOTE: Each core must have enough memory to store at least network-sized square matrix and given node sets to perform calculations.')	
	parser.add_argument('-bg', '--background', type=str, default='network', choices=['genesets', 'network'], required=False,
		help='Establishes the background gene set to calculate AUPRC over. Default is to use all genes in the network, can change to use only genes from the union of all gene sets tested (i.e. disease genes only).')	
	parser.add_argument('-ns', '--min_size', type=positive_int, default=20, required=False,
		help='Sets the minimum number of nodes in each gene set that must be found in the network to be evaluated. All node sets below this size are eliminated.')	
	parser.add_argument('-xs', '--max_size', type=positive_int, default=500, required=False,
		help='Sets the maximum number of nodes in each gene set that are be found in the network to be evaluated. All node sets above this size are eliminated.')	

	# Network performance score calculations (with null networks)
	parser.add_argument("-i", "--null_iter", type=positive_int, default=30, required=False,
		help='Number of times to perform degree-preserved shuffling of network to construct performance value null distribution. Default is 30. If this value is >0, --null_AUPRCs_save_path will be required')
	parser.add_argument('-nno', '--null_network_outdir', type=valid_outfile, default=None, required=False,
		help='File directory to save null networks after generation.')
	parser.add_argument('-nsp', '--null_AUPRCs_save_path', type=valid_outfile, default=None, required=False,
		help='CSV file path of where to save null network evaluation results. Used in the calculation of network performance score and perfomance gain scores')
	parser.add_argument('-psp', '--performance_save_path', type=valid_outfile, default=None, required=False,
		help='CSV file path of where to save network evaluation results as z-scores.')

	args = parser.parse_args()
	# If null networks need to be constructed
	if args.null_iter > 0:
		# A file path must be given to either save the null networks or the null network performance
		if (args.null_AUPRCs_save_path is None) and (args.null_network_outdir is None):
			parser.error('Save path required for null network edge lists or null network evaluation results.')

	####################################
	##### Network Evaluation Setup #####
	####################################

	# Limit core usage (if defined)
	import mkl
	mkl.set_num_threads(args.cores)
	
	# Load Network
	network = dit.load_network_file(args.network_path, verbose=args.verbose)
	network_size = len(network.nodes())

	# Load Gene sets
	genesets_raw = dit.load_node_sets(args.node_sets_file, verbose=args.verbose)

	# Calculate network kernel (also determine propagation constant if not set)
	kernel = nef.construct_prop_kernel(network, alpha=args.alpha, verbose=True)

	# Change background gene list if needed
	if args.background == 'genesets':
		background_node_set = set()
		for geneset in genesets:
			background_node_set = background_node_set.union(genesets[geneset])
		background_nodes = list(background_node_set.intersection(set(kernel.index)))
	else:
		background_nodes = list(kernel.index)

	# Filter gene sets
	print('Total node sets loaded:', len(genesets_raw))
	genesets = {gs:genesets_raw[gs].intersection(set(kernel.index)) for gs in genesets_raw 
	            if (len(genesets_raw[gs].intersection(set(kernel.index))) >= args.min_size) &
	               (len(genesets_raw[gs].intersection(set(kernel.index))) <= args.max_size)}
	print('Filtered node sets remaining:', len(genesets))

	# Calculate gene set sub-sample rate with network (if not set)
	if args.sample_p is None:
		genesets_p = nef.calculate_p(network, genesets)
	else:
		genesets_p = {geneset:args.sample_p for geneset in genesets}
	if args.verbose:
		print('Gene set sub-sample rates set')

	############################################
	##### Network Performance Calculations #####
	############################################

	# Calculate AUPRC for each gene set on actual network (large networks are >=10k nodes)
	if network_size < 10000:
		actual_AUPRC_values = nef.small_network_AUPRC_wrapper(kernel, genesets, genesets_p, n=args.sub_sample_iter, cores=args.cores, bg=background_nodes, verbose=False)
	else:
		actual_AUPRC_values = nef.large_network_AUPRC_wrapper(kernel, genesets, genesets_p, n=args.sub_sample_iter, cores=args.cores, bg=background_nodes, verbose=False)

	# Save the actual network's AUPRC values
	if args.verbose:
		print('Actual null network gene set AUPRCs calculated')
	actual_AUPRC_values.to_csv(args.actual_AUPRCs_save_path, header=False)

	#################################################
	##### Null Network Performance Calculations #####
	#################################################

	# If number of null networks > 0:
	if args.null_iter > 0:
		null_AUPRCs = []
		for i in range(args.null_iter):
			# Construct null networks and calculate AUPRCs for each gene set on each null network
			shuffNet = nef.shuffle_network(network, max_tries_n=10, verbose=True)
			# Save null network if null network output directory is given
			if args.null_network_outdir is not None:
				shuffNet_edges = shuffNet.edges()
				gct.write_edgelist(shuffNet_edges, args.null_network_outdir+'shuffNet_'+repr(i+1)+'.txt',
					delimiter='\t', binary=True)
				if args.verbose:
					print('Shuffled Network', i+1, 'written to file')
			# Construct null network kernel
			shuffNet_kernel = nef.construct_prop_kernel(shuffNet, alpha=args.alpha, verbose=False)
			# Calculate null network AUPRCs
			if network_size < 10000:
				shuffNet_AUPRCs = nef.small_network_AUPRC_wrapper(shuffNet_kernel, genesets, genesets_p, n=args.sub_sample_iter, cores=args.cores, bg=background_nodes, verbose=False)
			else:
				shuffNet_AUPRCs = nef.large_network_AUPRC_wrapper(shuffNet_kernel, genesets, genesets_p, n=args.sub_sample_iter, cores=args.cores, bg=background_nodes, verbose=False)
			null_AUPRCs.append(shuffNet_AUPRCs)
		# Construct table of null AUPRCs
		null_AUPRCs_table = pd.concat(null_AUPRCs, axis=1)
		null_AUPRCs_table.columns = ['shuffNet'+repr(i+1) for i in range(len(null_AUPRCs))]
		if args.verbose:
			print('All null network gene set AUPRCs calculated')
		# Save null network AUPRCs if save path is given
		if args.null_AUPRCs_save_path is not None:
			null_AUPRCs_table.to_csv(args.null_AUPRCs_save_path)
		# Calculate performance score/gain for each gene set's AUPRC if performance save path is given
		if args.performance_save_path is not None:
			network_gs_overlap = pd.Series({gs:len(genesets[gs]) for gs in genesets}, name='Genes in Network')
			network_performance = nef.calculate_network_performance_score(actual_AUPRC_values, null_AUPRCs_table, verbose=args.verbose)			
			network_perf_gain = nef.calculate_network_performance_gain(actual_AUPRC_values, null_AUPRCs_table, verbose=args.verbose)

			network_performance_table = pd.concat([network_gs_overlap, network_performance.rename('Performance Score'),
				                                   network_perf_gain.rename('Performance Gain')], axis=1, sort=False)
			network_performance_table.to_csv(args.performance_save_path)
		if args.verbose:
			print('Network performance score table saved')
