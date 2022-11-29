-- View: mvwchemmap_mapdrugbank

-- DROP MATERIALIZED VIEW mvwchemmap_mapdrugbank;

CREATE MATERIALIZED VIEW mvwchemmap_mapdrugbank
TABLESPACE pg_default
AS
 SELECT chem.drugbank_id,
    chem.smiles_clean,
    chem.inchikey,
    coo.dim1d2d,
    coo.dim3d,
    coo.neighbors_dim3,
    coo.neighbors_dimn,
    coo.desc_opera AS prop_value,
    prop.prop_value AS prop_tox,
    coo.d3_cube
   FROM chemical_description coo
     LEFT JOIN chemicals chem ON coo.inchikey = chem.inchikey
     LEFT JOIN chem_prop_drugbank_value prop ON chem.drugbank_id = prop.drugbank_id
  WHERE coo.map_name = 'drugbank'::text AND chem.drugbank_id IS NOT NULL
  GROUP BY chem.drugbank_id, chem.smiles_clean, chem.inchikey, coo.dim1d2d, coo.dim3d, coo.neighbors_dim3, coo.neighbors_dimn, prop.prop_value, coo.desc_opera, coo.d3_cube
WITH DATA;

ALTER TABLE mvwchemmap_mapdrugbank
    OWNER TO postgres;