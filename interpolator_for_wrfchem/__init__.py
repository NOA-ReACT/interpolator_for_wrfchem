import argparse

def main():
    """Entrypoint, parse arguments and run script"""

    argparser = argparse.ArgumentParser(description='Interpolate WRF-Chem output to a regular grid')
    argparser.add_argument('met_em_dir', help='Directory containing the met_em files')
    argparser.add_argument('input_file', help='Global model fields to interpolate')
    argparser.add_argument('wrfinput', help='WRF input file to update')
    argparser.add_argument('wrfbdy', help='WRF boundary file to update')
    args = argparser.parse_args()


