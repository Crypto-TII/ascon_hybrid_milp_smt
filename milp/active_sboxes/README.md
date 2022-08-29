# MILP model to compute the minimum number of active Sboxes for a given number of rounds.

## To run:
- run `make`
- `./a.out -t 12 -r 5`

Arguments:
- `-t`: number of threads
- `-r`: number of rounds

## Notes

- The current run will give a 5-round differential trail with 72 active Sboxes as an output.
- Modify the function "ascon_model" to verify other results.
- The file "input[3]_min[d2].txt gives the minimum number of active Sboxes at round 2 when
  there are only three active Sboxes at round 0 at positions (0, x, y).
