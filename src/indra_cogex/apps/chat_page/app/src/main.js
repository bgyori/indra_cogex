/* eslint-disable vue/multi-word-component-names */
import "./assets/indralab-style.css";
import "../../../static/source_badges.css"; // Generated after first import from indra_cogex/apps/constants.py

import { createApp, reactive } from "vue";

import { Statement } from "indralab-vue";
import { Evidence } from "indralab-vue";
import { CurationRow } from "indralab-vue";
import { RefLink } from "indralab-vue";
import { SourceDisplay } from "indralab-vue";
import App from "./App.vue";

const app = createApp(App);

app.component("Statement", Statement);
app.component("Evidence", Evidence);
app.component("RefLink", RefLink);
app.component("SourceDisplay", SourceDisplay);
app.component("CurationRow", CurationRow);

const cogexUrl = "https://discovery.indra.bio";

app.config.globalProperties.$stmt_hash_url = `${cogexUrl}/expand/`; // URL to get more evidence for a statement: e.g. https://discovery.indra.bio/expand/-16093215807632509?format=json-js&with_english=true&with_cur_counts=true&filter_ev=true
app.config.globalProperties.$curation_url = `${cogexUrl}/curate/`;
app.config.globalProperties.$curation_list_url = `${cogexUrl}/curation/list`;
app.config.globalProperties.$info_endpoint = `${cogexUrl}/chat/pusher_info`;
app.config.globalProperties.$sources = {
  readers: [
    "geneways",
    "tees",
    "gnbr",
    "isi",
    "trips",
    "rlimsp",
    "medscan",
    "eidos",
    "sparser",
    "reach",
  ],
  databases: [
    "psp",
    "cbn",
    "pc",
    "bel_lc",
    "signor",
    "biogrid",
    "tas",
    "hprd",
    "trrust",
    "ctd",
    "vhn",
    "pe",
    "drugbank",
    "omnipath",
    "conib",
    "crog",
    "dgi",
    "minerva",
    "creeds",
    "ubibrowser",
    "acsn",
  ],
};

const GStore = reactive({
  xrefs: {},
});

app.provide("GStore", GStore);
app.mount("#app");
