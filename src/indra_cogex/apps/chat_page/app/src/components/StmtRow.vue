<template>
  <div
    class="card-header stmt-row-header"
    :class="shown ? '' : 'border-bottom-0'"
    @click="shown = !shown"
    title="Click to show/hide statement details"
    data-toggle="collapse"
    :data-target="'#' + collapseAreaID"
    aria-expanded="false"
    :aria-controls="collapseAreaID"
  >
    <span>{{ stmt.english }}</span>
  </div>
  <div class="collapse" :id="collapseAreaID">
    <div class="card-body">
      <p class="text-center">Testing Statement</p>
      <Statement v-bind="fakeStatement" />
    </div>
  </div>
</template>

<script>
export default {
  name: "StmtRow.vue",
  data() {
    return {
      shown: false,
    };
  },
  props: {
    stmt: {
      type: Object,
      required: true,
    },
    rowIndex: {
      type: Number,
      required: true,
    },
  },
  computed: {
    fakeEvidence() {
      const fe = {
        text: "This is fake evidence text",
        pmid:
          this.stmt.evidence[0].pmid || this.stmt.evidence[0].text_refs.PMID,
        source_api: "sparser",
        text_refs: this.stmt.evidence[0].text_refs,
        num_curations: 0,
        source_hash: String(this.stmt.evidence[0].source_hash),
        stmt_hash: this.stmt.matches_hash,
        original_json: this.stmt.evidence[0],
      };
      console.log("fakeEvidence");
      console.log(fe);
      return fe;
    },
    fakeStatement() {
      const fs = {
        belief: 0.86,
        evidence: [this.fakeEvidence],
        english: this.stmt.english,
        hash: this.stmt.matches_hash,
        sources: { sparser: 1 },
        total_evidence: this.stmt.evidence.length,
        loadable: true,
      };
      console.log("fakeStatement");
      console.log(fs);
      return fs;
    },
    componentUUID() {
      return this.stmt.id;
    },
    collapseAreaID() {
      return `collapse-area-${this.stmt.id}`;
    },
    collapseControlID() {
      return `collapse-control-${this.stmt.id}`;
    },
  },
  methods: {
    getAgentList(stmt) {
      // Initialize agent list
      const list = [];
      // Loop through key-value pairs of the statement object and add
      // to agentList the objects that have the keys "name" and "db_refs"
      for (const obj of Object.values(stmt)) {
        if (obj.name && obj.db_refs) {
          list.push(obj);
        }
      }
      return list;
    },
  },
};
</script>

<style scoped>
.stmt-row-header {
  cursor: pointer;
}
</style>
