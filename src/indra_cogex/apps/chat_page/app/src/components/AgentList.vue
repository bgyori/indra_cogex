<template>
  <div class="text-left">
    Available agent types:
    <span
      v-for="([cls, descr], index) in availableClasses"
      class="badge"
      :title="availableClasses.length > 1 ? 'Click to toggle visibility' : ''"
      @click="toggleHide(cls)"
      :class="`${cls} ${isClsVisible(cls) ? '' : ' custom-opacity'} ${
        availableClasses.length > 1 ? 'hoverable' : ''
      }`"
      :key="index"
      >{{ descr }}</span
    >
  </div>
  <hr />
  <div class="text-left">
    <template v-for="(agent, number) in computedList" :key="number">
      <AgentModal :agent-object="agent" :button-class="'mx-1'" />
    </template>
  </div>
</template>

<script>
import AgentModal from "@/components/AgentModal.vue";

export default {
  name: "AgentList.vue",
  components: {
    AgentModal,
  },
  data() {
    return {
      badgeClasses: [
        ["badge-primary", "Gene/Protein"],
        ["badge-secondary", "Small Molecule"],
        ["badge-success", "Biological Process, Disease"],
        ["badge-info text-dark", "Phenotypic Abnormality"],
        ["badge-light text-dark", "Experimental Factor"],
        ["badge-warning text-dark", "ungrounded"],
      ],
      primaryVis: true,
      secondaryVis: true,
      successVis: true,
      infoVis: true,
      lightVis: true,
      warningVis: true,
    };
  },
  props: {
    entities: {
      type: Array,
      required: true,
    },
    showTitle: {
      // Show title
      type: Boolean,
      default: true,
    },
  },
  computed: {
    computedList() {
      let list = [];
      for (let agent of this.entities) {
        if (this.isClsVisible(this.dbRefsToCls(agent.db_refs))) {
          list.push(agent);
        }
      }
      return list;
    },
    availableClasses() {
      let classes = [];
      for (const agent of this.entities) {
        const cls = this.dbRefsToCls(agent.db_refs);
        // Add class if it is in badgeClasses and is not already in classes
        if (
          this.badgeClasses.some((el) => el[0] === cls) &&
          !classes.some((el) => el[0] === cls)
        ) {
          const clsDesc = this.badgeClasses.find((el) => el[0] === cls);
          classes.push(clsDesc);
        }
      }
      return classes;
    },
  },
  methods: {
    toggleHide(cls) {
      // If there is only one class, don't toggle
      if (this.availableClasses.length <= 1) {
        return;
      }
      switch (cls) {
        case "badge-primary":
          this.primaryVis = !this.primaryVis;
          break;
        case "badge-secondary":
          this.secondaryVis = !this.secondaryVis;
          break;
        case "badge-success":
          this.successVis = !this.successVis;
          break;
        case "badge-info text-dark":
          this.infoVis = !this.infoVis;
          break;
        case "badge-light text-dark":
          this.lightVis = !this.lightVis;
          break;
        case "badge-warning text-dark":
          this.warningVis = !this.warningVis;
          break;
      }
    },
    isClsVisible(cls) {
      switch (cls) {
        case "badge-primary":
          return this.primaryVis;
        case "badge-secondary":
          return this.secondaryVis;
        case "badge-success":
          return this.successVis;
        case "badge-info text-dark":
          return this.infoVis;
        case "badge-light text-dark":
          return this.lightVis;
        case "badge-warning text-dark":
          return this.warningVis;
      }
    },
    dbRefsToCls(db_refs) {
      // db_refs is an object with key-value pairs of the form
      // { "namespace": "namespace_id" }
      let cls = "";
      const dbRefKeys = Object.keys(db_refs);
      dbRefKeys.forEach((ns) => {
        // Continue as long as the class is not set or is still default
        if (cls === "warning text-dark" || cls === "") {
          const nsLower = ns.toLowerCase();
          switch (nsLower) {
            case "fplx":
            case "hgnc":
            case "up":
            case "uppro":
            case "mirbase":
              cls = "badge-primary";
              return;
            case "chebi":
              cls = "badge-secondary";
              return;
            case "go":
            case "mesh":
            case "doid":
              cls = "badge-success";
              return;
            case "hp":
              cls = "badge-info text-dark";
              return;
            case "efo":
              cls = "badge-dark";
              return;
            default:
              cls = "warning text-dark";
          }
        }
      });
      return cls;
    },
  },
};
</script>

<style scoped>
.hoverable {
  cursor: pointer;
}
.custom-opacity {
  opacity: 0.5;
}
</style>
