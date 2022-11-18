-- Table: public.tox21_assays

-- DROP TABLE public.tox21_assays;

CREATE TABLE public.tox21_assays
(
    protocol_name text COLLATE pg_catalog."default" NOT NULL,
    assay_target text COLLATE pg_catalog."default",
    cell_line text COLLATE pg_catalog."default",
    cell_target text COLLATE pg_catalog."default",
    description text COLLATE pg_catalog."default",
    design text COLLATE pg_catalog."default",
    CONSTRAINT protocol_name UNIQUE (protocol_name)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.tox21_assays
    OWNER to postgres;