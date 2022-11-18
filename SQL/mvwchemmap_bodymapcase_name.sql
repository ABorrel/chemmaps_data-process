-- View: public.mvwchemmap_bodymapcase_name

-- DROP MATERIALIZED VIEW public.mvwchemmap_bodymapcase_name;

CREATE MATERIALIZED VIEW public.mvwchemmap_bodymapcase_name
TABLESPACE pg_default
AS
 SELECT chemicals.casn,
    chemicals.name
   FROM chemicals
  WHERE (chemicals.casn IN ( SELECT bodymap_assay_ac50.casn
           FROM bodymap_assay_ac50))
WITH DATA;

ALTER TABLE public.mvwchemmap_bodymapcase_name
    OWNER TO postgres;