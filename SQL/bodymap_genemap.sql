-- Table: public.bodymap_genemap

-- DROP TABLE public.bodymap_genemap;

CREATE TABLE public.bodymap_genemap
(
    gene text COLLATE pg_catalog."default",
    organ text COLLATE pg_catalog."default",
    expression numeric,
    control numeric,
    system text COLLATE pg_catalog."default",
    id integer NOT NULL DEFAULT nextval('bodymap_genemap_1_id_seq'::regclass)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.bodymap_genemap
    OWNER to postgres;