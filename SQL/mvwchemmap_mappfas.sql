-- View: public.mvwchemmap_mappfas

-- DROP MATERIALIZED VIEW public.mvwchemmap_mappfas;

CREATE MATERIALIZED VIEW public.mvwchemmap_mappfas
TABLESPACE pg_default
AS
 SELECT chem.dsstox_id,
    chem.smiles_clean,
    chem.inchikey,
    coo.dim1d2d,
    coo.dim3d,
    coo.neighbors_dim3,
    coo.neighbors_dimn,
    coo.desc_opera AS prop_value,
    toxp.prop_value AS prop_tox,
    coo.d3_cube
   FROM chemical_description coo
     LEFT JOIN chemicals chem ON coo.inchikey = chem.inchikey
     LEFT JOIN chem_toxexp_value toxp ON toxp.dsstox_id = chem.dsstox_id
  WHERE coo.map_name = 'pfas'::text AND chem.dsstox_id IS NOT NULL AND coo.d3_cube IS NOT NULL
  GROUP BY chem.dsstox_id, chem.smiles_clean, chem.inchikey, coo.dim1d2d, coo.dim3d, coo.neighbors_dim3, coo.neighbors_dimn, coo.desc_opera, coo.d3_cube, toxp.prop_value
WITH DATA;

ALTER TABLE public.mvwchemmap_mappfas
    OWNER TO postgres;