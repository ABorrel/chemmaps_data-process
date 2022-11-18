-- Table: public.chemical_description

-- DROP TABLE public.chemical_description;
create sequence chemmap_3d_arr_id_seq;
CREATE TABLE public.chemical_description
(
    id bigint NOT NULL DEFAULT nextval('chemmap_3d_arr_id_seq'::regclass),
    source_id text COLLATE pg_catalog."default",
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
    CONSTRAINT chemmap_3d_arr_pkey PRIMARY KEY (id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.chemical_description
    OWNER to postgres;

COMMENT ON TABLE public.chemical_description
    IS 'Table includes dimension in 1d,2d and 3d also pre-calculated neighbors in 3d or nd. Can be searched by inchikey and source id (drugbank id ,dsstox-id)';

COMMENT ON COLUMN public.chemical_description.d3_cube
    IS 'set d3_cube[1] = dim1d2d[1],dim1d2d[2],dim3d[1]
Use to calculate Euclidean distance
';