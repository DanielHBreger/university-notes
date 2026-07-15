'''
Q4 entry point. The solution is split by exam section:
  q4_common.py - model functions and the shared fitting machinery
  q4_part_a.py - fits to the representative measurements
  q4_part_b.py - noise characterization of the real measurement
Each part file can also be run on its own.
'''
from q4_part_a import part_a
from q4_part_b import part_b
from q4_part_c import part_c

def main():
    file_path = 'Exam1.xlsx'
    part_a(file_path)
    part_b(file_path)
    part_c(file_path)


if __name__ == "__main__":
    main()
