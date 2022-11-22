-- Table: chem_descriptor_3d_name

-- DROP TABLE chem_descriptor_3d_name;

CREATE SEQUENCE chemmap_3d_arr_name_id_seq;
CREATE TABLE chem_descriptor_3d_name
(
    id bigint NOT NULL DEFAULT nextval('chemmap_3d_arr_name_id_seq'::regclass),
    name text COLLATE pg_catalog."default" NOT NULL,
    CONSTRAINT chemmap_3d_arr_name_pkey PRIMARY KEY (id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE chem_descriptor_3d_name
    OWNER to postgres;

COMMENT ON TABLE chem_descriptor_3d_name
    IS 'Table content the fields description for 3d array.  the order of the name column reflect the 3d array data stored in the chemmap_coords 
table column :dim3d';

GRANT USAGE, SELECT ON SEQUENCE chemmap_3d_arr_name_id_seq TO postgres;