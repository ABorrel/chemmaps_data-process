-- Table: chem_descriptor_opera_name

-- DROP TABLE chem_descriptor_opera_name;
create sequence chem_descriptor_opera_name_id_seq;
CREATE TABLE chem_descriptor_opera_name
(
    id integer NOT NULL DEFAULT nextval('chem_descriptor_opera_name_id_seq'::regclass),
    name text COLLATE pg_catalog."default",
    CONSTRAINT chem_descriptor_opera_name_pkey PRIMARY KEY (id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE chem_descriptor_opera_name
    OWNER to postgres;


GRANT USAGE, SELECT ON SEQUENCE chem_descriptor_opera_name_id_seq TO postgres;