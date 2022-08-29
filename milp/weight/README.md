# MILP model to compute the minimum weight of a r-round trail with a given configuration of Sboxes.

## To run:
- `make`
- `./a.out -t 12 -r 3`

## Notes
- The current run will produce a 3-round trail with weight 40 for the configuration [1, 3, 11].
- Change the configuration values in "int ID[3] = {1, 3, 11}" to find minimum weights for other 3-round configurations.

## Other files
- `W2.txt`: partial DDT table with weight 2. For instance the row (0, 0, 1, 0, 0, 0, 0, 1, 1, 0,, 1) means
   the input difference (0, 0, 1, 0, 0) goes to output difference (0, 0, 1, 1, 0) with probability $2^{-2}$ and
   has weight 2.
- `W3.txt`: partial DDT table with weight 3.
- `W4.txt`: partial DDT table with weight 4.
- `minimized_POS.txt`: The minimized product-of-sums (POS) formulas of above partial DDTs.
- `script_gen_ineq.py`: Generates vectors from minimized POS which are then used as a linear inequalities.
   To run: `python script_gen_ineq.py`
