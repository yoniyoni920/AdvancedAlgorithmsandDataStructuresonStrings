import sys
import abelian
import parametized
constants = []


def read_file(filename):
    try:
        with open(filename, 'r') as f:
            lines = [line.strip().upper() for line in f.readlines() if line.strip()]

        if "[" in lines[0]:
            num_sequences = int(lines[0].split()[-1])
        else:
            num_sequences = int(lines[0])

        return num_sequences, int(lines[1]), lines[2:]

    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)


def main():

    files = ['test1.txt', 'test2.txt', 'test3.txt']
    with open('output.txt', 'w', encoding='utf-8') as out_file:
        # Helper function to print and write at the same time
        def log(message):
            print(message)
            out_file.write(message + "\n")

        for filename in files:
            countmatches = 0
            num_lines, motif_length, text = read_file(filename)
            if not text:
                print("Error: Missing text data.")
                continue

            log(f"\n{'=' * 20}  {filename} {'=' * 20}")

            log("  Abelian  ")


            for k in range(num_lines):
                unique_patterns = set()
                flag = True
                for i in range(len(text[k]) - motif_length + 1):
                    pattern_p = text[k][i:i + motif_length]
                    unique_patterns.add(pattern_p)


                for pattern_p in unique_patterns:

                    result_obj = abelian.alebian(text[k], pattern_p)

                    if result_obj and result_obj['data']['matches']:

                        if flag:
                            log(f"\n[text] {result_obj['text']}")
                            flag = False
                        log(f"Found: {result_obj['data']['count']} pattern: '{pattern_p}':")


                        matches = result_obj['data']['matches']
                        for pos in matches:
                            found_text = text[k][pos: pos + motif_length]
                            log(f"  - Position {pos + 1}: Found '{found_text}'")
                        countmatches += len(matches)

            log(f" Found alebian'{countmatches}' matches in the text ")
            log(" ")
            log('=' * 100)
            log(" ")
            log(" Parameterized ")
            param_result_obj = parametized.discover_parameterized_results(text, motif_length, constants)

            if param_result_obj:

                 log(f"\n[text] {param_result_obj['text']}")

                 data = param_result_obj['data']
                 log(f"Structure Signature: {data['signature']}")
                 log(f"Common motif: '{data['motif']}'")
                 log(f"Total Occurrences: {data['count']}")
                 log("\nAll Results:")

                 for idx, match in enumerate(data['results']):
                    log(f"  Occurrence {idx + 1}:")
                    log(f"    Sequence: {match['sequence_num']}")
                    log(f"    Position: {match['position']}")
                    log(f"    Text: '{match['text']}'")
                    log(f"    Mapping (Swaps): {match['mapping']}")
                    log("-" * 30)
            else:
                 log("No Parameterized patterns found.")



if __name__ == "__main__":
    main()

