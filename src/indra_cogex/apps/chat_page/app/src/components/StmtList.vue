<template>
  <template v-if="metaIsLoaded">
    Available sources:
    <SourceDisplay :source-defaults="availableSources" />
  </template>
  <template v-else>
    <span>loading...</span>
    <div>
      <div class="spinner-grow spinner-grow-sm text-secondary" role="status">
        <span class="visually-hidden">Loading...</span>
      </div>
    </div>
  </template>
  <div v-if="!loading" class="list-group list-group-flush">
    <Statement
      v-for="(stmt, index) in stmtsWithEnglish"
      class="list-group-item"
      :key="`${this.componentUUID}-${index}`"
      :belief="stmt.belief"
      :evidence="stmt.evidence"
      :english="stmt.english"
      :hash="stmt.matches_hash"
      :sources="getSourceCounts(stmt.matches_hash)"
      :total_evidence="getEvidenceCount(stmt.matches_hash)"
      :loadable="true"
    />
  </div>
</template>

<script>
import helperFunctions from "@/helpers/helperFunctions";

export default {
  name: "StmtList.vue",
  data() {
    return {
      englishLookup: {},
      meta: {},
      loading: true,
    };
  },
  props: {
    stmts: {
      /* Array of statements:
       * Each statement is an object with the following properties:
       *  - belief: Number
       *  - evidence: Array of evidence objects {
       *      - annotations: Object containing evidence annotations
       *      - epistemics: Object containing evidence epistemics
       *      - pmid: String of article PMID
       *      - source_api: String of source API name (e.g. "sparser", "reach", "biopax", "signor")
       *      - source_hash: Number of source hash
       *      - text: String of evidence text
       *      - text_refs: Object containing identifiers for different article sources, e.g. { "PMID": "12345", "doi": "10.1038/ncomms12456" }
       *    }
       *  - id: Hex String UUID (what is this?)
       *  - matches_hash: String of statement hash
       *  - obj: (optional) Object containing the statement object as an Agent: { name: "", db_refs: {...} }
       *  - obj_activity: (optional) String describing the activity of the object
       *  - subj: (optional) Object containing the statement subject as an Agent: { name: "", db_refs: {...} }
       *  - subj_activity: (optional) String describing the activity of the subject
       *  - type: String describing the statement type (e.g. "Activation", "Inhibition")
       */
      type: Array,
      required: true,
    },
  },
  computed: {
    componentUUID() {
      return `stmtList-${this.uuid}`;
    },
    stmtsWithEnglish() {
      // Extend each statement object with its English assembled name
      if (Object.values(this.englishLookup).length === 0) {
        return this.stmts;
      }
      let extendedStmts = [];
      // Look up each statement in the English lookup
      this.stmts.forEach((stmt) => {
        const english = this.englishLookup[stmt.id];
        if (english) {
          stmt.english = english;
          extendedStmts.push(stmt);
        } else {
          extendedStmts.push(this.englishFromStmt(stmt));
        }
      });
      return extendedStmts;
    },
    availableSources() {
      let availableSources = { readers: [], databases: [] };

      // Loop the meta dict's values
      for (let metaEntry of Object.values(this.meta)) {
        // Loop the meta entry's keys of the sourceCounts dict
        for (let source of Object.keys(metaEntry.sourceCounts)) {
          // If the source is a reader, add it to the readers array
          const lowerSource = source.toLowerCase();

          if (
            this.$sources.readers.includes(lowerSource) &&
            !availableSources.readers.includes(lowerSource)
          ) {
            availableSources.readers.push(lowerSource);
          } else if (
            this.$sources.databases.includes(lowerSource) &&
            !availableSources.databases.includes(lowerSource)
          ) {
            availableSources.databases.push(lowerSource);
          } else if (
            !this.$sources.databases.includes(lowerSource) &&
            !this.$sources.readers.includes(lowerSource)
          ) {
            console.warn(`Unknown source: ${source}`);
          }
        }
      }
      return availableSources;
    },
    metaIsLoaded() {
      // Return true if the meta dict is not empty
      return Object.keys(this.meta).length > 0;
    },
  },
  mounted() {
    // Fetch English
    this.getEnglishArray();
    // Fetch statement metadata
    this.getStmtMeta();
    this.loading = false;
  },
  setup() {
    const uuid = helperFunctions.getUUID()();
    return {
      uuid,
    };
  },
  methods: {
    getAgentList(stmt) {
      // Initialize agent list
      const agentList = [];
      // Loop through key-value pairs of the statement object and add
      // to agentList the objects that have the keys "name" and "db_refs"
      for (const obj of Object.values(stmt)) {
        if (obj.name && obj.db_refs) {
          agentList.push(obj);
        }
      }
      return agentList;
    },
    stmtVerb(stmt) {
      // Return the verb of the statement noun, i.e. "activates" for "Activation"
      switch (stmt.type.toLowerCase()) {
        case "complex":
          return "binds";
        default:
          if (stmt.type.endsWith("ition")) {
            // If suffix is 'ition', replace with 'its'. E.g. 'Inhibition' -> 'Inhibits'
            return stmt.type.replace(/ition$/, "its");
          } else if (stmt.type.endsWith("ation")) {
            // If suffix is ation, replace with ates. E.g. "Activation" -> "Activates"
            return stmt.type.replace(/ation$/, "ates");
          } else {
            // If no suffix, return the type as-is
            return stmt.type.toLowerCase();
          }
      }
    },
    englishFromStmt(stmt) {
      // Backup in case there is no sentence in the English lookup
      const verb = this.stmtVerb(stmt);
      if (stmt.type.toLowerCase() === "complex") {
        // <first agent name> binds ', '.join(rest of agent names)
        const agentList = this.getAgentList(stmt);
        const firstAgent = agentList[0];
        const restOfAgents = agentList.slice(1);
        const restOfAgentsString = restOfAgents
          .map((agent) => agent.name)
          .join(", ");
        return `${firstAgent.name} ${verb} ${restOfAgentsString}`;
      }
      // <first agent name> <type> ', '.join(rest of agent names)
      const agentList = this.getAgentList(stmt);
      const firstAgent = agentList[0];
      const restOfAgents = agentList.slice(1);
      const restOfAgentsString = restOfAgents
        .map((agent) => agent.name)
        .join(", ");
      return `${firstAgent.name} ${verb} ${restOfAgentsString}`;
    },
    async getEnglishArray() {
      const indraApi = "http://api.indra.bio:8000/assemblers/english";
      // POST to get the English assembly of the statements
      const response = await fetch(indraApi, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          statements: this.stmts,
        }),
      });
      const data = await response.json();
      this.englishLookup = data.sentences;
    },
    getMeta(stmtHash) {
      if (this.meta[stmtHash]) {
        return this.meta[stmtHash];
      }
      return {};
    },
    getSourceCounts(stmtHash) {
      const meta = this.getMeta(stmtHash);
      if (meta.sourceCounts) {
        return meta.sourceCounts;
      }
      return {};
    },
    getEvidenceCount(stmtHash) {
      const meta = this.getMeta(stmtHash);
      if (meta.evidenceCount) {
        return meta.evidenceCount;
      }
      return 0;
    },
    getStmtMeta() {
      // Fetch statement metadata to https://discovery.indra.bio/api/get_stmts_meta_for_stmt_hashes
      const cogexApi =
        "https://discovery.indra.bio/api/get_stmts_meta_for_stmt_hashes";
      const stmtHashes = this.stmts.map((stmt) => stmt.matches_hash);
      fetch(cogexApi, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({
          stmt_hashes: stmtHashes,
        }),
      })
        .then((response) => response.json())
        .then((data) => {
          // Data should be list of relations in JSON format, map them to source counts and evidence counts per statement hash
          data.forEach((relation) => {
            const stmtObj = JSON.parse(relation.data.stmt_json);
            const stmtHash = stmtObj.matches_hash;
            const sourceCounts = JSON.parse(relation.data.source_counts);
            const evidenceCount = relation.data.evidence_count;
            this.meta[stmtHash] = {
              sourceCounts,
              evidenceCount,
            };
          });
        });
    },
  },
};
</script>
