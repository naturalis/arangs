import sys

input_file = sys.argv[1]

fasta_tab = open(input_file).readlines()

header = fasta_tab[0].strip()
sequence = ''.join(fasta_tab[1:]).replace('\n','')
start_stripped = sequence.lstrip('n')
end_stripped = sequence.rstrip('n')
cleaned_sequence = start_stripped.rstrip('n')
beginning_length = len(sequence) - len(start_stripped)
end_length = len(sequence) - len(end_stripped)
new_header = header + ' removed from beginning: ' + str(beginning_length)+ ', removed from end: ' + str(end_length)

print new_header
print cleaned_sequence

