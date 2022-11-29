-- Table: chemical_map

-- DROP TABLE chemical_map;
CREATE TABLE chemical_map
(
    inchikey text COLLATE pg_catalog."default",
    dim1d2d numeric[],
    dim3d numeric[],
    map_name text COLLATE pg_catalog."default",
    neighbors_dim3 text[] COLLATE pg_catalog."default",
    d3_cube numeric[],
    CONSTRAINT chem_on_map UNIQUE (inchikey,map_name)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE chemical_map
    OWNER to postgres;

COMMENT ON TABLE chemical_map
    IS 'Table includes chemical coordinates for each map in chemmaps';

COMMENT ON COLUMN chemical_map.d3_cube
    IS 'set d3_cube[1] = dim1d2d[1],dim1d2d[2],dim3d[1]
Use to calculate Euclidean distance
';