-- Table: chem_descriptor_1d2d_name

-- DROP TABLE chem_descriptor_1d2d_name;

create sequence chemmap_1d2d_arr_name_id_seq;

CREATE TABLE chem_descriptor_1d2d_name
(
    id bigint NOT NULL DEFAULT nextval('chemmap_1d2d_arr_name_id_seq'::regclass),
    name text COLLATE pg_catalog."default" NOT NULL,
    CONSTRAINT chemmap_1d2d_arr_name_pkey PRIMARY KEY (id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE chem_descriptor_1d2d_name
    OWNER to postgres;

COMMENT ON TABLE chem_descriptor_1d2d_name
    IS 'Table content the fields description for 2d array.  the order of the name column matches the 1d2d array data stored in the chemmap_coords 
table dim1d2d ';

GRANT USAGE, SELECT ON SEQUENCE chemmap_1d2d_arr_name_id_seq TO postgres;
