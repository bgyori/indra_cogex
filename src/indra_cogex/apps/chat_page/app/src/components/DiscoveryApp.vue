<template>
  <h3>Discovery App</h3>
  <!-- Chat app -->
  <template v-if="logged_in">
    <p>Logged in as {{ chat.name }} ({{ chat.email }})</p>
    <div class="input-group mb-3">
      <input
        type="text"
        class="form-control"
        placeholder="Enter question"
        aria-label="Enter question"
        aria-describedby="button-addon"
        v-model="text_input"
        :disabled="disable_input"
        @keyup.enter="sendMessage"
      />
      <button
        class="btn btn-outline-secondary"
        type="button"
        id="button-addon"
        @click="sendMessage"
        :disabled="disable_input"
      >
        Ask
      </button>
    </div>
    <div class="clearfix message d-flex flex-column-reverse">
      <template v-for="(message, index) in chat.history" :key="index">
        <MessageWrapper
          :user="message.user"
          :bot="message.bot"
          @clarification-requested2="(text) => (this.text_input = text)"
        />
      </template>
    </div>
  </template>
  <!-- Login -->
  <template v-else>
    <form id="loginForm" @submit.prevent="logIntoChatSession">
      <div class="form-row">
        <div class="col-md-6">
          <div class="form-group">
            <label for="name" class="form-label">Name</label>
            <input
              type="text"
              class="form-control"
              id="name"
              placeholder="Name*"
              v-model="form_name"
              :disabled="submitted"
              required
            />
          </div>
        </div>
        <div class="col-md-6">
          <div class="form-group">
            <label for="email" class="form-label">Email</label>
            <input
              type="email"
              class="form-control"
              id="email"
              placeholder="Email*"
              v-model="form_email"
              :disabled="submitted"
              required
            />
          </div>
        </div>
      </div>
      <button type="submit" class="btn btn-primary">Start Session</button>
    </form>
  </template>
</template>

<script>
import Pusher from "pusher-js";
import MessageWrapper from "@/components/MessageWrapper.vue";

export default {
  name: "DiscoveryApp.vue",
  components: { MessageWrapper },
  props: {
    info_endpoint: {
      type: String,
      default: "/api_chat/pusher_info", // Set corresponding endpoint in indra_cogex/apps/chat_page/chat_blueprint.py
    },
  },
  data() {
    return {
      /** Pusher options and related data
       * {
       *    "pusher_key": pusher_key or "",
       *    "pusher_cluster": pusher_cluster or "",
       *    "auth_endpoint": url_for(".pusher_authentication"),
       *    "new_user_endpoint": url_for(".guestUser"),
       * }
       */
      pusher_info: {},
      pusher: null,
      form_name: "",
      form_email: "",
      chat: {
        name: undefined,
        email: undefined,
        channel: undefined,
        response: undefined, // Last response from chatbot
        history: [], // Store query and response: [{user: {...}, bot: {...}}, ...]
      },
      submitted: false,
      logged_in: false,
      disable_input: false,
      text_input: "",
    };
  },
  mounted() {
    this.getAppInfo();
  },
  computed: {
    pusher_key() {
      if (this.pusher_info) {
        return this.pusher_info.pusher_key;
      }
      return "";
    },
    pusher_cluster() {
      if (this.pusher_info) {
        return this.pusher_info.pusher_cluster;
      }
      return "";
    },
    auth_endpoint() {
      if (this.pusher_info && this.pusher_info.auth_endpoint) {
        return this.pusher_info.auth_endpoint;
      }
      return "";
    },
    new_user_endpoint() {
      if (this.pusher_info && this.pusher_info.new_user_endpoint) {
        return this.pusher_info.new_user_endpoint;
      }
      return "";
    },
  },
  methods: {
    async getAppInfo() {
      const resp = await fetch(this.info_endpoint, {
        method: "GET",
      });
      const data = await resp.json();
      this.pusher_info = await data;
      // Setup pusher instance
      this.pusher = await this.setupPusher();
      // Bind to message events
      this.pusher.bind("client-support-new-message", (data) => {
        this.newMessage(data);
      });
      // Bind to chunked message events
      let events = {};
      this.pusher.bind("chunked-client-support-new-message", (data) => {
        if (!Object.prototype.hasOwnProperty.call(events, data.id)) {
          events[data.id] = { chunks: [], receivedFinal: false };
        }
        let ev = events[data.id];
        ev.chunks[data.index] = data.chunk;
        if (data.final) ev.receivedFinal = true;
        if (
          ev.receivedFinal &&
          ev.chunks.length === Object.keys(ev.chunks).length
        ) {
          console.log("Received chunked message");
          let msgJson = JSON.parse(ev.chunks.join(""));
          this.newMessage(msgJson);
          delete events[data.id];
        }
      });
    },
    async setupPusher() {
      return new Pusher(this.pusher_key, {
        authEndpoint: this.auth_endpoint,
        cluster: this.pusher_cluster,
        encrypted: true,
      });
    },
    async logIntoChatSession() {
      if (this.form_email && this.form_name) {
        // Disable the form
        this.submitted = true;
        // Trim strings
        let trimmed_name = this.form_name.trim();
        let trimmed_email = this.form_email.trim();
        const resp = await fetch(this.new_user_endpoint, {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
          },
          body: JSON.stringify({
            name: trimmed_name,
            email: trimmed_email,
          }),
        });
        const data = await resp.json();
        const user_info = await data;
        // Set the chat session info
        this.chat.name = trimmed_name;
        this.chat.email = trimmed_email;
        // Subscribe to the channel
        this.chat.channel = this.pusher.subscribe("private-" + user_info.email);
        this.logged_in = true;
      }
    },
    sendMessage() {
      if (this.chat.channel) {
        // Disable the input
        this.disable_input = true;

        // Create the message
        const createdAt = new Date().toUTCString();
        const message = {
          sender: this.chat.name,
          email: this.chat.email,
          text: this.text_input,
          createdAt: createdAt,
        };

        // Send the message to the channel
        this.chat.channel.trigger("client-guest-new-message", message);

        // Add the message to the message history
        let user_msg = {
          input: message.text,
          name: message.sender,
          sender_email: message.email,
          createdAt: message.createdAt,
          responseIndex: undefined, // Set index to undefined until response is received
        };

        // Push message to history
        this.chat.history.push({ user: user_msg, bot: { receivedAt: null } });

        // Clear the input
        this.text_input = "";

        // Enable the input
        this.disable_input = false;
        console.log("Message sent");
      } else {
        console.log("No channel, cannot send message!");
      }
    },
    async newMessage(message) {
      // Expecting at least {raw_text: "", name: "", objects: ""}
      console.log("New message received; bot=");
      console.log(message);
      if (message !== undefined) {
        // await the message
        let resolved_message = await message;
        let receivedAt = new Date().toUTCString();
        // Replace the response to the message history
        this.chat.history[this.chat.history.length - 1].bot = {
          ...resolved_message,
          receivedAt: receivedAt,
        };
        console.log("Last interaction after response push: ");
        console.log(this.chat.history[this.chat.history.length - 1]);
      }
    },
  },
};
</script>
