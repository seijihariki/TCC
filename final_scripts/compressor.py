#!/usr/bin/python3

import math as m

__DEBUG__ = False


def compress(input_it):
    'Compresses input lines to output lines (stream compressor)'

    def output(start_value, end_value, seq_length):
        out = f"{start_value}"
        if end_value and end_value != start_value:
            out += f"-{end_value}"
        if seq_length and seq_length != 1:
            out += f"x{seq_length}"
        return out

    # Current and last two number streaks (value, length, part of streak)
    number_streaks = [[None, None, False]] * 3

    # Sequence data
    sequence_start_value = None
    sequence_start_direction = None

    # Extra executions for finishing properly
    extra = 2

    input_it = iter(input_it)

    while extra > 0:
        value = None
        try:
            value = int(next(input_it))
        except StopIteration:
            extra -= 1

        # If a new number streak has started
        if not value or value != number_streaks[0][0]:
            if __DEBUG__:
                print(number_streaks)

            # Finalize sequence if unable to continue, like when streak
            # changes size or step size
            if sequence_start_value and (
                    number_streaks[0][1] != number_streaks[1][1] or
                    sequence_start_direction != number_streaks[0][0] - number_streaks[1][0]):
                out = output(sequence_start_value,
                             number_streaks[1][0], number_streaks[1][1])
                if __DEBUG__:
                    print(f'Writing {out}')
                yield out
                sequence_start_value = None
                sequence_start_direction = None
                number_streaks[1][2] = True
            # Start sequence if found
            elif not sequence_start_value and\
                    number_streaks[0][1] and\
                    number_streaks[0][0] and\
                    number_streaks[0][1] == number_streaks[1][1] and\
                    abs(number_streaks[0][0] - number_streaks[1][0]) == 1:
                sequence_start_value = number_streaks[1][0]
                sequence_start_direction = number_streaks[0][0] - \
                    number_streaks[1][0]
                number_streaks[1][2] = True
                if __DEBUG__:
                    print(
                        f"Starting sequence from {sequence_start_value}[{sequence_start_direction}]")
            # Set as sequence
            elif sequence_start_value:
                number_streaks[1][2] = True

            # If unique repeated value (Not part of a sequence)
            if not number_streaks[1][2] and number_streaks[1][0]:
                out = output(number_streaks[1][0],
                             None, number_streaks[1][1])
                if __DEBUG__:
                    print(f'Writing {out}u')
                yield out

            if __DEBUG__:
                print(f"Starting {value} streak")
            # Update number streaks
            number_streaks[2] = number_streaks[1]
            number_streaks[1] = number_streaks[0]
            number_streaks[0] = [value, 0, False]

        if value:
            number_streaks[0][1] += 1


def decompress(input_it):
    'Decompresses input lines to output lines (stream compressor)'

    it = iter(input_it)

    for line in it:
        # Parse compressed file

        # Split line
        range_values, *seq_length_v = line.split('x')
        start_s, *end_v = range_values.split('-')

        # Parse values
        start = int(start_s)
        end = int(end_v[0]) if len(end_v) == 1 else start
        seq_length = int(seq_length_v[0]) if len(seq_length_v) == 1 else 1

        # Decompress
        step = 1 if end > start else - 1
        for v in range(abs(end - start) + 1):
            current_val = start + v * step
            for _ in range(seq_length):
                yield current_val


# For command line use
if __name__ == "__main__":
    import argparse

    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Compresses and Decompresses number sequence files")

    parser.add_argument('filename', help='Input Filename')

    parser.add_argument("-d", "--decompress",
                        help="Decompress file", action="store_true")

    parser.add_argument("-o", "--output-file",
                        type=str, help="Output File")

    parser.add_argument("-v", "--verbose",
                        help="Print debug Information", action="store_true")

    args = parser.parse_args()

    __DEBUG__ = args.verbose

    if args.decompress:
        if not args.filename.endswith(".cseq"):
            print(f"File must end with .cseq")

        output_file = args.output_file or args.filename[-5]

        print(f"Decompressing file: {args.filename} into {output_file}")
        with open(output_file, 'w') as output:
            output.writelines(
                map(lambda x: f"{x}\n", decompress(open(args.filename))))
    else:
        output_file = args.output_file or f"{args.filename}.cseq"

        if not output_file.endswith(".cseq"):
            print(f"File must end with .cseq")

        print(f"Compressing file: {args.filename} into {output_file}")
        with open(output_file, 'w') as output:
            output.writelines(
                map(lambda x: f"{x}\n", compress(open(args.filename))))
