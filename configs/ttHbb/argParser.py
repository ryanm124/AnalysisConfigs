import argparse
import os,getpass

def get_year_from_args():
    parser = argparse.ArgumentParser(description='Run analysis on NanoAOD files using PocketCoffea processors')
    # Inputs
    parser.add_argument("-c","--cfg", default=os.getcwd() + "/config/test.py", required=True, type=str,
                        help='Config file with parameters specific to the current run')
    parser.add_argument("-ro", "--custom-run-options", type=str, default=None, help="User provided run options .yaml file")
    parser.add_argument("-o", "--outputdir", required=True, type=str, help="Output folder")
    parser.add_argument("-y", "--year", required=True, type=str, help="year")
    parser.add_argument("-sa", "--sample", required=True, type=str, nargs='+', help="sample name")
    parser.add_argument("-t", "--test", action="store_true", help="Run with limit 1 interactively")
    parser.add_argument("-lf","--limit-files", type=int, help="Limit number of files")
    parser.add_argument("-lc","--limit-chunks", type=int, help="Limit number of chunks", default=None)
    parser.add_argument("-e","--executor", type=str,
                        help="Overwrite executor from config (to be used only with the --test options)" )
    parser.add_argument("-s","--scaleout", type=int, help="Overwrite scalout config" )
    parser.add_argument("-ll","--loglevel", type=str, help="Logging level", default="INFO" )
    parser.add_argument("-f","--full", action="store_true", help="Process all datasets at the same time", default=False )
    args = parser.parse_args()
    return [args.year,args.sample ]
