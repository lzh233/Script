import argparse
from script.log import add_log
from script.snp_sumarrize import Sumarrize_Variant

@add_log
def main():
    parser = argparse.ArgumentParser(description='Get the summarize report of a vid')
    parser.add_argument('--data_dir', help='The directory which save the results of snp', required=True)
    parser.add_argument('--results_dir', help='The results directory to save the vid sumarrize', required=True)
    parser.add_argument('--detect_share', action="store_true",help='Find the vid which in both alt and ref', required=False)
    parser.add_argument('--manual_sample_list', action="store_true",help='Manual make a sample list', required=False)
    parser.add_argument('--sample_list_dir',help='location of manual sample list', required=False)

    args = parser.parse_args()
    run = Sumarrize_Variant(data_dir = str(args.data_dir),
                            results_dir=str(args.results_dir),
                            detect_share = args.detect_share,
                            manual_sample_list = args.manual_sample_list,
                            sample_list_dir = str(args.sample_list_dir))
    
    run.get_vid_summarize()
    #os.system("rm -rf sample_list.txt")


if __name__ == '__main__':
    main()
