from readfasta import readfasta
import glob, os

def main():
    print( "*****\nBioinformatics - Assignment 3 - Group 2\n*****\n" )

    # scan in all fasta files in the "genes" directory
    os.chdir( os.getcwd() + "/genes/" )
    file = glob.glob( "*.fasta" )[0]

    # read all the genes from the fasta file
    print( file )
    genes = readfasta( file )
    print( genes )

main()

