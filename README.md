# Compute map database
Scripts used to compute map coordinate and manage database used for ChemMaps
- 08-08-2019: Create README
- 08-09-2019: Update

            -Update folder split folder in R and py
            -Convert in python3.6 format (0%)
            -Add descriptor py files from server
            -Convert in python3.6 (50%)

- 08-15-2019: finish convertion python
- 09-04-2019: Update for the new DrugBank and remove descriptors scripts, used the library developed to generalize the descriptor computation
- 09-05-2019: Clean source for drugbank and dell sources for descriptor computation, part of the generalisation protocol
- 09-06-2019: Add SQL class to interact with the database
- 09-12-2019: Add SQL sequence in the drugbank class and optimize run on the DSSTOX
- 09-20-2019: Update DSSTOX lib for pfas and Tox21
- 09-24-2019: connect to DB for transfer
- 09-26-2019: compute the map split ans centroid, optimize fonction split and add function to load only first coords
- 10-08-2019: Add update map on the dsstox for the split


# Todo list
- 09-20-2019: finish to upload prop table in the DSSTOX table
- 09-26-2019: compute the neighbor for DSSTOX -> done for n = 3
- 09-26-2019: update png file
