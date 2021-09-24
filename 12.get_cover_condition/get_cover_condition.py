import argparse
from script.log import add_log
from script.statistics_bam import Static_Bam

@add_log
def main():

    parser = argparse.ArgumentParser(description='Get the summarize report of a vid')
    parser.add_argument('--data_dir', help='The directory which save the data', required=True)
    parser.add_argument('--results_dir', help='The results directory to save the vid sumarrize', required=True)
    parser.add_argument('--manual_sample_list', action="store_true",help='Manual make a sample list', required=False)
    parser.add_argument('--sample_list',help='Location of manual sample list', required=False)
    parser.add_argument('--bp_min',help = "Provide a min length to filter data(default: 60)", required=False,default = 50)

    args = parser.parse_args()
    run = Static_Bam(data_dir = str(args.data_dir),
                            results_dir=str(args.results_dir),
                            manual_sample_list = args.manual_sample_list,
                            sample_list_dir = str(args.sample_list),
                            bp_min = int(args.bp_min))
    
    run.get_cover_condition()
    
if __name__ == '__main__':
    main()
