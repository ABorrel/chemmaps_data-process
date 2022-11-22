-- Table: chemical_description

-- DROP TABLE chemical_description;
CREATE TABLE chemical_description
(
    inchikey text COLLATE pg_catalog."default",
    dim1d2d numeric[],
    dim3d numeric[],
    map_name text COLLATE pg_catalog."default",
    neighbors_dim3 text[] COLLATE pg_catalog."default",
    neighbors_dimn text[] COLLATE pg_catalog."default",
    d3_cube numeric[],
    desc_1d2d numeric[],
    desc_3d numeric[],
    interference_prediction numeric[],
    desc_opera numeric[],
    CONSTRAINT chem_in_map UNIQUE (inchikey,map_name)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE chemical_description
    OWNER to postgres;

COMMENT ON TABLE chemical_description
    IS 'Table includes dimension in 1d,2d and 3d also pre-calculated neighbors in 3d or nd. Can be searched by inchikey and source id (drugbank id ,dsstox-id)';

COMMENT ON COLUMN chemical_description.d3_cube
    IS 'set d3_cube[1] = dim1d2d[1],dim1d2d[2],dim3d[1]
Use to calculate Euclidean distance
';