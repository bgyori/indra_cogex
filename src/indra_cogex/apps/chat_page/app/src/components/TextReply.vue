<template>
  <!-- Show spinner if receivedDate is null -->
  <div v-if="!receivedDate" class="d-flex align-items-center">
    <span>Loading...</span>
    <div
      class="spinner-grow spinner-grow-sm text-secondary ml-auto"
      role="status"
      aria-hidden="true"
    ></div>
  </div>
  <template v-if="raw_text">
    <div class="text-left">
      <template
        v-for="(txtObj, index) in textObjectArray"
        :key="`${componentID}-text${index}`"
      >
        <!-- No object: just display the text as innerHtml -->
        <span v-if="txtObj.object === null" v-html="txtObj.text"></span>
        <!-- A single agent from query: display as AgentModal -->
        <template
          v-else-if="
            txtObj.object.role === 'fixed' && txtObj.object.type === 'agent'
          "
        >
          <AgentModal :agent-object="txtObj.object.value" />
        </template>
        <!-- A string list, set as an unordered list if  -->
        <template v-else-if="txtObj.object.type === 'string_list'">
          <template v-if="txtObj.format === 'bullet'">
            <ul>
              <li
                v-for="(str, index) in txtObj.object.value"
                :key="`${componentID}-strlist${index}`"
              >
                {{ str }}
              </li>
            </ul>
          </template>
          <!-- Assumed to be a pre-formatted flat list -->
          <template v-else>
            <span v-html="txtObj.text"></span>
          </template>
        </template>
        <!-- agent list with reply entities: display the computed string -->
        <span
          v-else-if="
            txtObj.object.type === 'agent_list' &&
            txtObj.object.role === 'variable'
          "
          v-html="txtObj.text"
        ></span>
        <!-- URL list -->
        <span v-else-if="txtObj.object.type === 'url_list'">
          <span
            v-for="([url, db_name], index) in txtObj.object.value"
            :key="`${componentID}-urllist${index}`"
          >
            <a :href="url" target="_blank">{{ db_name }}</a>
            <span v-if="index < txtObj.object.value.length - 2">, </span>
            <span v-else-if="index === txtObj.object.value.length - 2"
              >, and
            </span>
          </span>
        </span>
        <!-- Clarification: display the computed string -->
        <div v-if="clarification">
          <button
            class="btn btn-sm btn-outline-secondary"
            title="Click to copy the clarification text to the input text box"
            @click="$emit('clarification-requested1', clarification)"
          >
            Use clarification
          </button>
        </div>
      </template>
    </div>
  </template>
</template>

<script>
import AgentModal from "@/components/AgentModal.vue";
import helperFunctions from "@/helpers/helperFunctions";

export default {
  name: "TextReply.vue",
  emits: ["clarification-requested1"],
  components: {
    AgentModal,
  },
  props: {
    raw_text: {
      type: String,
      default: "",
    },
    objects: {
      type: Object,
      default: null,
    },
    receivedDate: {
      type: String,
      default: "",
    },
  },
  computed: {
    clarification() {
      // Match 'Your question is similar to "%s". Try asking it that way.'
      const clarificationRegex =
        /Your question is similar to "(.*)". Try asking it that way./;
      const match = clarificationRegex.exec(this.raw_text);
      if (match) {
        return match[1];
      }
      return "";
    },
    componentID() {
      return `text-reply-${this.uuid}`;
    },
    textObjectArray() {
      // Locate all references to text objects in raw text

      // Contains objects with entries: "text", "object"
      // "text" is the text to show, "object" is the object representing the text (null if plain text)
      let txtObjs = [];
      const rawText = this.raw_text || "";

      // Split string on any text between '{' and '}'
      let txtObjStrings = rawText.split(/(\{.*?\})/g);

      if (txtObjStrings.length > 0 && this.objects !== null) {
        // Loop through each string and match against objects
        txtObjStrings.forEach((substring) => {
          // Skip empty strings
          if (substring.length > 0) {
            // If we match a recognised object, look it up among the objects
            if (substring.includes("{") && substring.includes("}")) {
              let objName = substring.replace(/[{}]/g, "");
              let formatSpec;
              // Check if there is a format specifier
              if (objName.includes(":")) {
                [objName, formatSpec] = objName.split(":");
              } else {
                formatSpec = null;
              }
              // Get the actual object
              const obj = this.objects[objName] || null;
              let textForObj;
              if (obj) {
                textForObj = this.getTextForObj(obj, formatSpec);
                txtObjs.push({
                  text: textForObj,
                  object: obj,
                  objName: objName,
                  match: substring, // The string that matched the object
                  format: formatSpec, // The format specifier
                });
              } else {
                console.warn(
                  `No object available for object name ${objName}`,
                  this.objects
                );
              }
            } else {
              // If we don't match a recognised object, just add the string
              txtObjs.push({
                text: substring,
                object: null,
                objName: null,
                match: substring,
                format: null,
              });
            }
          }
        });
      } else {
        if (rawText.includes("{") || rawText.includes("}")) {
          console.warn(
            "MessageWrapper.vue: No text objects found in raw text, but references found. Returning raw text.",
            rawText
          );
        }
        txtObjs.push({
          text: this.raw_text,
          object: null,
          objName: null,
          match: null,
          format: null,
        });
      }
      return txtObjs;
    },
  },
  methods: {
    englishJoin(arr, delimiter = ", ", oxford = true) {
      // Join an array of strings into a string with a delimiter and a grammatically correct conjunction
      if (arr.length === 0) {
        console.warn("No array to join", arr);
      } else if (arr.length === 1) {
        return arr[0];
      } else if (arr.length === 2) {
        return `${arr[0]} and ${arr[1]}`;
      } else {
        return `${arr.slice(0, -1).join(delimiter)}${
          oxford ? "," : ""
        } and ${arr.slice(-1)}`;
      }
    },
    getFormatSpec(formatSpec) {
      if (!formatSpec) {
        return null;
      }
      // Get the format specifier from a string
      if (formatSpec.includes("limit")) {
        // Match limit(d+) and return the number
        const limitRegex = /limit(\d+)/;
        const match = limitRegex.exec(formatSpec);
        if (match) {
          return parseInt(match[1]);
        }
        console.warn("Couldn't match limit specifier", formatSpec);
        return null;
      } else {
        return formatSpec;
      }
    },
    getTextForObj(obj, formatSpec = null) {
      // Return the text for the object
      const formatSpecVal = this.getFormatSpec(formatSpec);
      switch (obj.type) {
        case "agent":
          // Return object.value.name
          return obj.value.name;
        case "agent_list": {
          // If there is a limit, return the first N agent names as a string
          let limit =
            formatSpec && formatSpec.includes("limit")
              ? formatSpecVal
              : obj.value.length;
          let nameArray = obj.value.slice(0, limit).map((agent) => agent.name);
          return this.englishJoin(nameArray);
        }
        case "string_list":
          /* string_list is an array of strings that are joined together either
           *  as a bullet list (handled in template) or by comma (done here w englishJoin) */
          if (formatSpecVal === "bullet") return null;
          return this.englishJoin(obj.value);

        case "url_list":
          return null; // Always show handled in template
        case "string":
          // Return the string
          return obj.value;
        default:
          console.log("In switch-case; unhandled object type: " + obj.type);
          console.log(obj);
          return "(N/A)";
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
