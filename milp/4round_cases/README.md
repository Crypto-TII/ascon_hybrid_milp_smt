# MILP model for the number of configurations which may give a 4-round trail with at most 44 active Sboxes.

## To run:
- `python script.py`

## Notes
- `possible_configurations.txt`: Full list of possible configurations. Obtained by the MILP model.
   To run the MILP model, simply type:

  - ` make`
  - `./a.out -t 12`

- `infeasible_3r.txt`: List of 3-round infeasible configurations.
- The output of `python script.py` is a list of reduced number of configurations.
