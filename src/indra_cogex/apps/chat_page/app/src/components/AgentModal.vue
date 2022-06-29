<template>
  <!-- BootStrap 4.6 Modal -->
  <!-- Button trigger modal -->
  <button
    type="button"
    :title="title"
    @click.once="fillXrefs()"
    class="btn badge"
    :class="`${badgeClass} ${buttonClass ? buttonClass : ''}`"
    data-toggle="modal"
    :data-target="`#${modalUUID}`"
  >
    {{ badgeText }}
  </button>

  <!-- Modal -->
  <div
    class="modal fade"
    :id="modalUUID"
    tabindex="-1"
    :aria-labelledby="titleUUID"
    aria-hidden="true"
  >
    <div class="modal-dialog">
      <div class="modal-content">
        <div class="modal-header">
          <h5 class="modal-title" :id="titleUUID">{{ textToShow }}</h5>
          <button
            type="button"
            class="close"
            data-dismiss="modal"
            aria-label="Close"
          >
            <span aria-hidden="true">&times;</span>
          </button>
        </div>
        <div class="modal-body">
          <p>
            <span v-if="entityDescription" v-html="entityDescription"></span
            ><span v-else><i>Description missing</i></span>
          </p>
          <hr />
          <table class="table table-striped">
            <thead>
              <tr>
                <th>Database</th>
                <th>Identifier</th>
                <th>Lookup</th>
              </tr>
            </thead>
            <tbody>
              <tr
                v-for="(xref, index) in allRefs"
                :key="`${modalUUID}-row${index}`"
              >
                <td>{{ xref[0] }}</td>
                <td>{{ xref[1] }}</td>
                <td>
                  <a :href="`${xref[2]}`" target="_blank">
                    <i class="fas fa-external-link-alt"></i>
                  </a>
                </td>
              </tr>
            </tbody>
          </table>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import { DefaultValues } from "@/helpers/DefaultValues";
import helperFunctions from "@/helpers/helperFunctions";

export default {
  name: "AgentModal.vue",
  inject: ["GStore"],
  props: {
    text: {
      type: String,
      default: null,
    },
    agentObject: {
      // Object of { name: "", db_refs: {ns: id} }
      type: Object,
      required: true,
      validator: (obj) => {
        return obj.name !== null && obj.db_refs !== null;
      },
    },
    buttonClass: {
      // Optional. Class to apply to the button for easier styling
      type: String,
      default: null,
    },
  },
  data() {
    return {
      xrefs: [],
      serverError: false,
      lookupData: {},
    };
  },
  computed: {
    badgeText() {
      let text = this.text || this.agentObject.name;
      const maxLength = 30;
      const half = Math.floor(maxLength / 2);
      if (text.length > maxLength) {
        text =
          text.substring(0, half) + "..." + text.substring(text.length - half);
      }
      return text;
    },
    textToShow() {
      return (
        this.text || (this.agentObject ? this.agentObject.name : "(no name)")
      );
    },
    topGrounding() {
      // Set to the first grounding in db_refs
      let [firstNs, firstId] = Object.entries(this.agentObject.db_refs)[0];
      let topEntry = [firstNs.toLowerCase(), firstId];
      // Return the top ranked grounding according to the nsPrio map, lower is better
      for (const [ns, id] of Object.entries(this.agentObject.db_refs)) {
        const nsLower = ns.toLowerCase();
        if (this.getNsPrio(nsLower) < this.getNsPrio(topEntry[0])) {
          topEntry = [nsLower, id];
        }
      }
      if (
        DefaultValues.nsPriorityMap[topEntry[0].toLowerCase()] === undefined
      ) {
        console.warn(
          `${topEntry[0]} is not in the nsPrio map. Forgot to make it lowercase in code?`
        );
      }
      return topEntry;
    },
    title() {
      return `Grounded to ${this.topGrounding.join(":")}. Click to see more.`;
    },
    modalUUID() {
      return `modal-${this.uuid}`;
    },
    titleUUID() {
      return `modal-${this.uuid}-title`;
    },
    allRefs() {
      if (this.topGrounding.length === 0) {
        return this.xrefs;
      }

      return [
        [
          this.topGrounding[0],
          this.topGrounding[1],
          `https://bioregistry.io/${this.topGrounding[0]}:${this.topGrounding[1]}`,
        ],
        ...this.xrefs,
      ];
    },
    entityDescription() {
      if (this.lookupData.definition) {
        return this.lookupData.definition
          .replace(/\\n/g, "<br>")
          .replace(/\n/g, "<br>");
      }
      return "";
    },
    badgeClass() {
      for (const [cls, values] of Object.entries(DefaultValues.badgeMappings)) {
        if (
          this.topGrounding[0] &&
          values.includes(this.topGrounding[0].toLowerCase())
        ) {
          return cls;
        }
      }
      console.log(
        `No badge class found for ${this.textToShow} (${this.topGrounding.join(
          ":"
        )})`
      );
      return "text-dark border";
    },
  },
  methods: {
    // todo: Add fetch of synonyms from http://grounding.indra.bio/get_names POST endpoint
    async fillXrefs() {
      // Call network search web api to get xrefs; fixme: should use a standalone api; Isn't there one for the bioontology?
      if (this.topGrounding.length === 0) {
        console.log(
          "Cannot fill crossrefs without grounding for",
          this.textToShow
        );
        return;
      }
      const topNsUpper = this.topGrounding[0].toUpperCase();
      if (this.xrefs.length <= 1) {
        const key = this.topGrounding.join(":");
        // Only fetch if we don't have any xrefs yet
        if (this.GStore.xrefs && this.GStore.xrefs[key]) {
          this.xrefs = this.GStore.xrefs[key];
          console.log(`Fetched xrefs for ${key} from cache`);
        } else {
          const xrefsUrl = `https://network.indra.bio/api/xrefs?ns=${topNsUpper}&id=${this.topGrounding[1]}`;
          const xrefResp = await fetch(xrefsUrl);
          const xrefData = await xrefResp.json();
          this.xrefs = await xrefData;
          this.GStore.xrefs[key] = this.xrefs;
        }
      } else {
        console.log(`Already have xrefs for ${this.textToShow}`);
      }

      if (Object.entries(this.lookupData).length === 0) {
        // Call biolookup.io, e.g. http://biolookup.io/api/lookup/DOID:14330
        const bioluUrl = `http://biolookup.io/api/lookup/${topNsUpper}:${this.topGrounding[1]}`; // Currently only supports http
        const bioluResp = await fetch(bioluUrl);
        const bioluData = await bioluResp.json();
        this.lookupData = await bioluData;
      }
    },
    getNsPrio(ns) {
      const nsPrio = DefaultValues.nsPriorityMap;
      if (!ns) {
        return nsPrio.default;
      } else if (nsPrio[ns] === undefined) {
        return nsPrio.default;
      } else {
        return nsPrio[ns];
      }
    },
  },
  setup() {
    const uuid = helperFunctions.getUUID()();
    return {
      uuid,
    };
  },
};
</script>
