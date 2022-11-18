-- Table: public.chem_descriptor_opera_name_new

-- DROP TABLE public.chem_descriptor_opera_name_new;

CREATE TABLE public.chem_descriptor_opera_name_new
(
    id integer NOT NULL,
    name text COLLATE pg_catalog."default",
    CONSTRAINT chem_descriptor_opera_name_new_pkey PRIMARY KEY (id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.chem_descriptor_opera_name_new
    OWNER to postgres;