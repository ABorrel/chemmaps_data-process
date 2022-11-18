-- Table: public.chemical_description_user

-- DROP TABLE public.chemical_description_user;

CREATE TABLE public.chemical_description_user
(
    id bigint NOT NULL DEFAULT nextval('chemmap_coord_user_id_seq'::regclass),
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
    status text COLLATE pg_catalog."default",
    CONSTRAINT id PRIMARY KEY (id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.chemical_description_user
    OWNER to postgres;