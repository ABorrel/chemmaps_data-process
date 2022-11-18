-- Table: public.chem_toxexp_name

-- DROP TABLE public.chem_toxexp_name;

CREATE TABLE public.chem_toxexp_name
(
    id bigint NOT NULL,
    name text COLLATE pg_catalog."default",
    CONSTRAINT chem_toxexp_name_pkey PRIMARY KEY (id)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.chem_toxexp_name
    OWNER to postgres;