-- Table: public.bodymap_assay_ac50

-- DROP TABLE public.bodymap_assay_ac50;

CREATE TABLE public.bodymap_assay_ac50
(
    casn text COLLATE pg_catalog."default",
    assay text COLLATE pg_catalog."default",
    ac50 numeric
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.bodymap_assay_ac50
    OWNER to postgres;