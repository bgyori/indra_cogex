<template>
  <div v-if="user && user.input" :id="componentID" class="row">
    <div class="col card mb-1 px-0">
      <div
        class="card-header msg-wrapper-header"
        :class="shown ? '' : 'border-bottom-0'"
        title="Click to expand/collapse response"
        @click="shown = !shown"
        data-toggle="collapse"
        :data-target="`#${idRegistry.collapseID}`"
        aria-expanded="true"
        :aria-controls="idRegistry.collapseID"
      >
        <div class="row">
          <div class="col-auto">
            <span class="text-muted small" :title="user.createdAt"
              >User Input</span
            >
          </div>
          <div class="col text-left">
            <span>{{ user.input }}</span>
          </div>
        </div>
      </div>
      <div class="collapse show" :id="idRegistry.collapseID">
        <div class="card-body">
          <div class="row mb-1">
            <div class="col-2 text-left">
              <span
                v-if="receivedDate"
                class="text-muted small"
                :title="receivedDate"
                >Output ({{ shortDate }})
              </span>
            </div>
            <div class="col-10"></div>
          </div>
          <div class="row">
            <!-- Awrapping col-12 is needed to not screw up the tabs -->
            <div class="col-12">
              <!-- tabs -->
              <nav>
                <div
                  class="nav nav-tabs"
                  :id="idRegistry.navTabsID"
                  role="tablist"
                >
                  <!-- Text tab -->
                  <a
                    v-if="bot && bot.raw_text"
                    class="nav-link active"
                    title="See text response from bot"
                    :id="idRegistry.nav.textID"
                    data-toggle="tab"
                    :href="`#${idRegistry.content.textID}`"
                    role="tab"
                    :aria-controls="idRegistry.content.textID"
                    aria-selected="true"
                    >Response</a
                  >
                  <!-- Entities tab -->
                  <a
                    v-if="replyEntities.length > 0"
                    class="nav-link"
                    @click.once="click.entities = true"
                    title="See agents associated with this response"
                    :id="idRegistry.nav.entitiesID"
                    data-toggle="tab"
                    :href="`#${idRegistry.content.entitiesID}`"
                    role="tab"
                    :aria-controls="idRegistry.content.entitiesID"
                    aria-selected="false"
                    >Agents ({{ replyEntities.length }})</a
                  >
                  <!-- Stmts tab -->
                  <a
                    v-if="replyStmts.length > 0"
                    class="nav-link"
                    @click.once="click.stmts = true"
                    title="See statements associated with this response"
                    :id="idRegistry.nav.stmtsID"
                    data-toggle="tab"
                    :href="`#${idRegistry.content.stmtsID}`"
                    role="tab"
                    :aria-controls="idRegistry.content.stmtsID"
                    aria-selected="false"
                    >Statements ({{ replyStmts.length }})</a
                  >
                </div>
              </nav>
              <!-- Tab content -->
              <div class="tab-content" :id="idRegistry.tabContentID">
                <!-- Text tab -->
                <div
                  class="tab-pane fade show active overflow-container"
                  :id="idRegistry.content.textID"
                  role="tabpanel"
                  :aria-labelledby="idRegistry.nav.textID"
                >
                  <div class="card card-body border-light">
                    <TextReply
                      :raw_text="bot ? bot.raw_text : ''"
                      :objects="objects"
                      :received-date="receivedDate"
                      @clarification-requested1="
                        this.emitClarificationRequested
                      "
                    />
                  </div>
                </div>
                <!-- Entities content -->
                <div
                  v-if="replyEntities.length > 0"
                  class="tab-pane fade overflow-container"
                  :id="idRegistry.content.entitiesID"
                  role="tabpanel"
                  :aria-labelledby="idRegistry.nav.entitiesID"
                >
                  <div class="card card-body border-light">
                    <template v-if="click.entities">
                      <AgentList :entities="replyEntities" />
                    </template>
                  </div>
                </div>
                <!-- Stmts content -->
                <div
                  v-if="replyStmts.length > 0"
                  class="tab-pane fade overflow-container"
                  :id="idRegistry.content.stmtsID"
                  role="tabpanel"
                  :aria-labelledby="idRegistry.nav.stmtsID"
                >
                  <div class="card card-body border-light">
                    <template v-if="click.stmts">
                      <StmtList :stmts="replyStmts" />
                    </template>
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import AgentList from "@/components/AgentList.vue";
import StmtList from "@/components/StmtList.vue";
import TextReply from "@/components/TextReply.vue";
import helperFunctions from "@/helpers/helperFunctions";

export default {
  name: "MessageWrapper.vue",
  emits: ["clarification-requested2"],
  components: {
    TextReply,
    AgentList,
    StmtList,
  },
  data() {
    return {
      listMax: 0,
      shown: true,
      click: {
        text: false,
        entities: false,
      },
    };
  },
  props: {
    bot: {
      /* Expecting a bot object:
       * {
       *  raw_text: "",
       *  name: "",
       *  objects: {...},
       *  sender: "",
       *  receivedAt: "",
       * }
       * */
      type: Object,
      default: null,
    },
    user: {
      /* Expecting this object:
       * {
       *  input: "",
       *  name: "",
       *  sender_email: "",
       *  createdAt: "",
       * }
       */
      type: Object,
      default: null,
    },
  },
  computed: {
    receivedDate() {
      if (this.bot.receivedAt !== null) {
        return new Date(this.bot.receivedAt).toLocaleString();
      }
      return "";
    },
    shortDate() {
      // Format: "HH:mm" 24 hour clock
      if (this.bot.receivedAt !== null) {
        return new Date(this.bot.receivedAt).toLocaleTimeString("en-GB", {
          hour: "2-digit",
          minute: "2-digit",
        });
      }
      return "";
    },
    objects() {
      if (this.bot && this.bot.objects) {
        return this.bot.objects;
      }
      return null;
    },
    replyEntities() {
      // Loop through the objects and check if they have role === 'variable'
      if (this.bot && this.bot.objects) {
        return this.getObjectsWithTypeRole(
          this.bot.objects,
          "agent_list",
          "variable"
        );
      }
      return [];
    },
    replyStmts() {
      // Loop through the objects and check if they have role === 'statement'
      if (this.bot && this.bot.objects) {
        return this.getObjectsWithTypeRole(
          this.bot.objects,
          "statement_list",
          "variable"
        );
      }
      return [];
    },
    componentID() {
      return `msgWrapper-${this.uuid}`;
    },
    idRegistry() {
      return {
        collapseID: `${this.componentID}-collapse`,
        navTabsID: `${this.componentID}-nav-tabs`,
        tabContentID: `${this.componentID}-tabs-content`,
        nav: {
          textID: `${this.componentID}-tab-text`,
          entitiesID: `${this.componentID}-tab-entities`,
          stmtsID: `${this.componentID}-tab-stmts`,
        },
        content: {
          textID: `${this.componentID}-content-text`,
          entitiesID: `${this.componentID}-content-entities`,
          stmtsID: `${this.componentID}-content-stmts`,
        },
      };
    },
  },
  methods: {
    isActive(tab) {
      // Check if the tab is active
      return this.activeTab === tab;
    },
    getObjectsWithTypeRole(objs, type, role) {
      // Get all entities from type 'agent_list' with a given role
      let entities = [];
      // Loop the entries in the Object
      Object.values(objs).forEach((obj) => {
        if (obj && obj.type === type && obj.role === role) {
          // Concatenate the entities (stored in the 'value' field)
          entities = entities.concat(obj.value);
        }
      });
      return entities;
    },
    emitClarificationRequested(text) {
      // Emit the event to the parent component
      this.$emit("clarification-requested2", text);
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

<style scoped>
.overflow-container {
  max-height: 300px;
  overflow-y: auto;
}
.msg-wrapper-header {
  cursor: pointer;
}
</style>
