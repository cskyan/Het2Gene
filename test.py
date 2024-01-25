import argparse
parser = argparse.ArgumentParser(description='Demo of argparse')
parser.add_argument('--contents', type=str, default="Het2GeneTest")

args=parser.parse_args()
contents = args.contents

print("helloword")

