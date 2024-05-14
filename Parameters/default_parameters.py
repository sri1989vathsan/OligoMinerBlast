#Allow user to input parameters on command line.
# userInput = argparse.ArgumentParser(description=\
#     '%s version %s. Requires a FASTA file as input. Currently, only '
#     'single-entry FASTA files are supported.  Returns a .fastq file, which '
#     'can be inputted into short read alignment programs. Optionally, a '
#     '.bed file can be outputted instead if \'-b\' is flagged. Tm values '
#     'are corrected for [Na+] and [formamide].' % (scriptName, Version))
# requiredNamed = userInput.add_argument_group('required arguments')
# requiredNamed.add_argument('-f', '--file', action='store', required=True,
#                             help='The FASTA file to find probes in')
# userInput.add_argument('-l', '--minLength', action='store', default=36,
#                         type=int,
#                         help='The minimum allowed probe length; default is '
#                             '36')
# userInput.add_argument('-L', '--maxLength', action='store', default=41,
#                         type=int,
#                         help='The maximum allowed probe length, default is '
#                             '41')
# userInput.add_argument('-g', '--min_GC', action='store', default=20,
#                         type=int,
#                         help='The minimum allowed percent G + C, default is '
#                             '20')
# userInput.add_argument('-G', '--max_GC', action='store', default=80,
#                         type=int,
#                         help='The maximum allowed percent  G + C, default '
#                             'is 80')
# userInput.add_argument('-t', '--min_Tm', action='store', default=42,
#                         type=int,
#                         help='The minimum allowed Tm, default is 42')
# userInput.add_argument('-T', '--max_Tm', action='store', default=47,
#                         type=int,
#                         help='The maximum allowed Tm, default is 47')
# userInput.add_argument('-X', '--prohibitedSeqs', action='store',
#                         default='AAAAA,TTTTT,CCCCC,GGGGG', type=str,
#                         help='Prohibited sequence list (separated by commas '
#                             'with no spaces), default is '
#                             '\'AAAAA,TTTTT,CCCCC,GGGGG\'')
# userInput.add_argument('-s', '--salt', action='store', default=390,
#                         type=int,
#                         help='The mM Na+ concentration, default is 390')
# userInput.add_argument('-F', '--formamide', action='store', default=50,
#                         type=float,
#                         help='The percent formamide being used, default is '
#                             '50')
# userInput.add_argument('-S', '--Spacing', action='store', default=0,
#                         type=int,
#                         help='The minimum spacing between adjacent probes, '
#                             'default is 0 bases')
# userInput.add_argument('-c', '--dnac1', action='store', default=25,
#                         type=float,
#                         help='Concentration of higher concentration strand '
#                             '[nM] -typically the probe- to use for '
#                             'thermodynamic calculations. Default is 25')
# userInput.add_argument('-C', '--dnac2', action='store', default=25,
#                         type=float,
#                         help='Concentration of lower concentration strand '
#                             '[nM] -typically the target- to use for '
#                             'thermodynamic calculations. Default is 25')
# userInput.add_argument('-n', '--nn_table', action='store',
#                         default='DNA_NN3',
#                         type=str,
#                         help='The nearest neighbor table of thermodynamic '
#                             'parameters to be used. See options in '
#                             'Bio.SeqUtils.MeltingTemp. Default is DNA_NN3')
# userInput.add_argument('-H', '--header', action='store', type=str,
#                         help='Allows the use of a custom header in the '
#                             'format chr:start-stop. E.g. '
#                             '\'chr2:12500-13500\'')
# userInput.add_argument('-b', '--bed', action='store_true', default=False,
#                         help='Output a .bed file of candidate probes '
#                             'instead of a .fastq file.')
# userInput.add_argument('-O', '--OverlapMode', action='store_true',
#                         default=False,
#                         help='Turn on Overlap Mode, which returns all '
#                             'possible candidate probes in a block of '
#                             'sequence including overlaps. Off by default. '
#                             'Note, if selecting this option, the '
#                             '-S/--Spacing value will be ignored')
# userInput.add_argument('-v', '--verbose', action='store_true',
#                         default=False,
#                         help='Turn on verbose mode to have probe mining'
#                             'progress print to Terminal. Off by default')
# userInput.add_argument('-R', '--Report', action='store_true', default=False,
#                         help='Write a Report file detailing the results of '
#                             'each window of sequence considered by the '
#                             'script. The first set of lines give the '
#                             'occurrence of each possible failure mode for '
#                             'quick reference. Off by default. Note, '
#                             'selecting this option will slow the script '
#                             'considerably')
# userInput.add_argument('-D', '--Debug', action='store_true', default=False,
#                         help='The same as -Report, but prints info to '
#                             'terminal instead of writing a log file. Off '
#                             'by default')
# userInput.add_argument('-M', '--Meta', action='store_true', default=False,
#                         help='Write a text file containing meta information '
#                             'Off by default. Reports input file <tab> '
#                             'estimated runtime <tab> blockParse version '
#                             '<tab> candidates discovered <tab> span in kb '
#                             'covered by candidate probes <tab> candidate '
#                             'probes per kb')
# userInput.add_argument('-o', '--output', action='store', default=None,
#                         type=str, help='Specify the stem of the output '
#                                         'filename')


# 
# Import user-specified command line values.
# args = userInput.parse_args()
# inputFile = args.file
l = 36
L = 41
gcPercent = 20
GCPercent = 80
tm = 42
TM = 47
X = 'AAAAA,TTTTT,CCCCC,GGGGG'
sal = 390
form = 50
sp = 0
concA = 25
concB = 25
headerVal = None
bedVal = False
OverlapModeVal = False
verbocity = False
reportVal = False
debugVal = False
metaVal = False
outNameVal = None
nntable='DNA_NN3'